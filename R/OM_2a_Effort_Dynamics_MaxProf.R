#-------------------------------------------------------------------------------
# Functions to be used within nloptr to obtain effort and effortshare according
# to maximization of profits, restringed to Fcube like approach in 
# TAC constraints.
# restrictions on effort by fleet:
#
# Dorleta Garcia and Agurtzane Urtizberea - Azti Tecnalia
# Created: 03/06/2016
# Changed: 23/11/2011 10:16:01 - Sonia Sanchez
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# MaxProfit(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1)
# Same as original MaxProfit +
#  + min and max effort thresholds by fleet
#  + info on discrimination capability of the metiers (fleets.ctrl[[fl]]$q2zero[st,mt])
#-------------------------------------------------------------------------------
MaxProfit <- function(fleets, biols, BDs,covars, advice, biols.ctrl, fleets.ctrl, advice.ctrl, flnm, year = 1, season = 1,...){
  
  dimnms <- dimnames(biols[[1]]@n)
  
  if(length(year) > 1 | length(season) > 1)
    stop('Only one year and season is allowed' )
  
  # If year/season/iter numerics => indicate position 
  # else names => get positions.
  if(length(year) > 1 | length(season) > 1)
    stop('Only one year and season is allowed' )
  
  # 'year' dimension.
  yr <- year
  if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
  if(length(yr) == 0) stop('The year is outside object time range')  
  
  # 'season' dimension.
  ss <- season
  if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
  if(length(ss) == 0) stop('The season is outside object season range')  
  
  
  # Dimensions.
  nst <- length(biols);          stnms <- names(biols)
  ns  <- dim(biols[[1]]@n)[4]
  it  <- dim(biols[[1]]@n)[6]
  flnms <- names(fleets)
  
  # Fleet's info
  fl    <- fleets[[flnm]]
  sts   <- intersect(fleets.ctrl[[flnm]][['stocks.restr']], catchNames(fl))
  mtnms <- names(fl@metiers)
  nmt   <- length(mtnms)
  
  Et.res  <-  numeric(it)
  names(Et.res) <- dimnms$iter
  efs.res <-  matrix(NA,nmt,it,dimnames = list(mtnms, dimnms$iter))
  
  
  # Check fleets.ctrl elements.
  # This checks the restriction for all the fleets and we only need to check the fleet we are proyecting.
  #   if(! all(sapply(names(fleets), function(x) fleets.ctrl[[x]]$restriction %in% c('catch', 'landings'))))
  #       stop("fleets.ctrl$restriction must be equal to 'catch' or 'landings'")
  
  effort.restr <- ifelse(length(fleets.ctrl[[flnm]]$effort.restr) > 1, fleets.ctrl[[flnm]]$effort.restr[yr], fleets.ctrl[[flnm]]$effort.restr)
  restriction <- ifelse(length(fleets.ctrl[[flnm]]$restriction) > 1,   fleets.ctrl[[flnm]]$restriction[yr], fleets.ctrl[[flnm]]$restriction)
  if(!(restriction %in% c('catch', 'landings') )) stop("fleets.ctrl$restriction for fleet, ', flnm, ', must be equal to 'catch' or 'landings'")
  
  LO <- ifelse(length(fleets.ctrl[[flnm]]$LandObl)>1, fleets.ctrl[[flnm]]$LandObl[yr], fleets.ctrl[[flnm]]$LandObl) 
  
  # The effort is restricted only by the stocks in 'stocks.restr'                      
  if(!is.null(fleets.ctrl[[flnm]][['stocks.restr']])) {
    stocks.restr <- fleets.ctrl[[flnm]][['stocks.restr']]
  } else
    stocks.restr <- fleets.ctrl[[flnm]][['stocks.restr']] <- catchNames(fleets[[flnm]])
  
  # if(is.na(fleets.ctrl[[flnm]][['stocks.restr']]))
  #   cat(paste("warning: fleets.ctrl[['",flnm,"']][['stocks.restr']]==NA, then effort is restricted by capacity and not by any of the stocks.\n", sep=""))
  
  
  for(i in dimnms$iter){
    
    # For avoiding errors in R CMD CHECK: 
    # add variables that will be defined within FLObjs2S3_fleetSTD function call
    B <- N <- QS <- TAC <- rho <- efs.m <- vc.m <- fc <- crewS <- effs <- Cr.f <- TAC.yr <- 
      tacos <- q.m <- alpha.m <- beta.m <- pr.m <- ret.m <- wd.m <- wl.m <- K <- 
      Nyr_1 <- Myr_1 <- M <- Cfyr_1 <- Cyr_1 <- LO <- NULL
    
    # Advice season for each stock
    adv.ss <- setNames( rep(NA,nst), stnms)
    for (st in stnms) adv.ss[st] <- ifelse(is.null(advice.ctrl[[st]][["adv.season"]]), ns, advice.ctrl[[st]][["adv.season"]]) # [nst]
    
    # Transform the FLR objects into list of arrays in order to be able to work with non-FLR
    list2env(FLObjs2S3_fleetSTD(biols = biols, fleets = fleets, advice = advice, covars = covars, 
                                biols.ctrl = biols.ctrl, fleets.ctrl = fleets.ctrl, BDs=BDs, 
                                flnm = flnm, yr = yr, ss = ss, iters = i, adv.ss = adv.ss), environment())
    
    
    # # Correction of efs.m when no TAC for main target species
    if(!is.null(fleets.ctrl[[flnm]]$q2zero)){ 
      fl.sel <- fleets.ctrl[[flnm]]$q2zero
      for (mt in mtnms)
        if (sum(Cr.f[rownames(fl.sel)[fl.sel[,mt]==0 & !is.na(fl.sel[,mt])],])==0) { 
          efs.m[mtnms!=mt,] <- efs.m[mtnms!=mt,] + efs.m[mtnms!=mt,] * efs.m[mt,]/sum(efs.m[mtnms!=mt,])
          efs.m[mt,] <- 0
        }
      if(fleets.ctrl[[flnm]]$efs.abs == FALSE) efs.m <- efs.m/sum(efs.m)
    }
    
    # Calculate the initial point based on the effort that correspond with the TAC quotas.
    effs <- numeric(length(q.m))
    names(effs) <- names(q.m)
    # Numbers at age by stock
    Nsts <- lapply(setNames(sts, sts), function(x) array(N[[x]][,,i,drop=T], dim = c(dim(N[[x]])[1:2],1)))
    
    for(st in names(q.m)){
      
      effort.fun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'effort', sep = '.')

      
      effs[st] <-  eval(call(effort.fun, Cr = Cr.f,  N = N, q.m = q.m, rho = rho, efs.m = efs.m, 
                                alpha.m = alpha.m, beta.m = beta.m, ret.m = ret.m, wl.m = wl.m, wd.m = wd.m,stknm=st,
                                restriction = restriction,  QS.groups = fleets.ctrl[[flnm]][['QS.groups']],
                                tac=TAC[st,i, drop=F], Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1,  M = M, Cfyr_1 = Cfyr_1))
      
    }
    
    qsum.stk <- sapply(names(q.m), function(x) sum(q.m[[x]]))
    
    Et  <- 0.9*ifelse(effort.restr == 'min', min(effs[qsum.stk>0]), effs[effort.restr])[1] 
    Et <- ifelse(Et < K, Et, K*0.9)
    
    
    catch.restr <- ifelse(is.null(restriction), 'landings', restriction)
    
    # if(is.null(fleets.ctrl[[flnm]]$opts)) opts <- list("algorithm" = "NLOPT_LN_COBYLA", maxeval = 1e9, xtol_abs = rep(1e-4,nmt), xtol_rel = 1e-4, maxtime = 300)
    # else  opts <- fleets.ctrl[[flnm]]$opts
    
    # Effort restrictions (by metier)
    if(!is.null(fleets.ctrl[[flnm]]$effort.range)){
      effort.range <- fleets.ctrl[[flnm]]$effort.range
      efs.max <- as.numeric(effort.range[,"max"])
      efs.min <- as.numeric(effort.range[,"min"])
    }
    else{
      if(fleets.ctrl[[flnm]]$efs.abs == FALSE){
        efs.min <- rep(0, length(fleets[[flnm]]@metiers))
        efs.max <- rep(1, length(fleets[[flnm]]@metiers))
        names(efs.min) <- names(efs.max) <- names(fleets[[flnm]]@metiers)
      }else{
        efs.min <- rep(0, length(fleets[[flnm]]@metiers))
        efs.max <- rep(K, length(fleets[[flnm]]@metiers))
        names(efs.min) <- names(efs.max) <- names(fleets[[flnm]]@metiers)
      }
    }
    
    #  if(fleets.ctrl[[flnm]]$efs.abs == FALSE){
    # If efs.min == 0 or Cr.f == 0 or E0 == 0 set them equal to 1e-8 to avoid having indeterminations in the penalties
    E0 <- Et*efs.m
    
    efs.min <- ifelse(efs.min == 0, 1e-8, efs.min)
    E0      <- ifelse(E0 == 0, 1e-8, E0)
    Cr.f[Cr.f == 0] <-  1e-8
    
    # recalculate efs.m in case E0 has changed.
    efs.m <- E0/sum(E0)
    
    # Apply these restrictions to initial values
    if (fleets.ctrl[[flnm]]$efs.abs == FALSE) {
      
      efs.min <- ifelse(efs.m <= efs.min, efs.m*0.99, efs.min)
      efs.m <- ifelse(efs.m >= efs.max, efs.max*0.99, efs.m)
      
      E0 <- ifelse(Et*efs.m == 0, 1e-8, Et*efs.m)
      
    } else {
      
      efs.min <- ifelse(E0 <= efs.min, E0*0.99, efs.min)
      E0 <- ifelse(E0 >= efs.max, efs.max*0.99, E0)
      
    }
    
    # }
    #  else{
    #    E0 <- sum(efs.m
    #  }
    
    X <- log(E0/(K - E0))
    
    L <- list(X=X,E0=E0, effs=effs, efs.max = efs.max, efs.min = efs.min,q.m = q.m, alpha.m = alpha.m, 
              beta.m = beta.m, pr.m = pr.m, ret.m = ret.m, wd.m = wd.m,
              wl.m = wl.m, N = N, B = B, fc = fc, vc.m = vc.m,   Cr.f = Cr.f,  crewS = crewS, K = K , 
              effort.restr = effort.restr, catch.restr = catch.restr, stocks.restr = stocks.restr, efs.abs = fleets.ctrl[[flnm]]$efs.abs, 
              tacos = tacos, rho = rho, tac=TAC[,i, drop= F], Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1,  M = M, Cfyr_1 = Cfyr_1)

    eff_opt <- try(optim(X, f_MP_nloptr_penalized, efs.max = efs.max, 
      efs.min = efs.min, q.m = q.m, alpha.m = alpha.m, 
      beta.m = beta.m, pr.m = pr.m, ret.m = ret.m, wd.m = wd.m, 
      wl.m = wl.m, N = N, B = B, fc = fc, vc.m = vc.m, 
      Cr.f = Cr.f, crewS = crewS, K = K, effort.restr = effort.restr, 
      catch.restr = catch.restr, stocks.restr = stocks.restr, 
      efs.abs = fleets.ctrl[[flnm]]$efs.abs, tacos = tacos, 
      rho = rho, tac = TAC[, i, drop = F], Cyr_1 = Cyr_1, 
      Nyr_1 = Nyr_1, Myr_1 = Myr_1, M = M, Cfyr_1 = Cfyr_1, 
      flnm = flnm, fleets.ctrl = fleets.ctrl), silent = TRUE)
    
    if(class(eff_opt) != "try-error"){ #! not recommended use in ?try --> Use if(inherits(res, "try-error")) instead
      if(eff_opt[["convergence"]] %in% c(1, 10)){
        eff_opt <- try(optim(eff_opt[["par"]], f_MP_nloptr_penalized, 
          efs.max = efs.max, efs.min = efs.min, q.m = q.m, 
          alpha.m = alpha.m, beta.m = beta.m, pr.m = pr.m, 
          ret.m = ret.m, wd.m = wd.m, wl.m = wl.m, N = N, 
          B = B, fc = fc, vc.m = vc.m, Cr.f = Cr.f, crewS = crewS, 
          K = K, effort.restr = effort.restr, catch.restr = catch.restr, 
          stocks.restr = stocks.restr, efs.abs = fleets.ctrl[[flnm]]$efs.abs, 
          tacos = tacos, rho = rho, tac = TAC[, i, drop = F], 
          Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1, 
          M = M, Cfyr_1 = Cfyr_1, flnm = flnm, fleets.ctrl = fleets.ctrl, 
          control = list(maxit = 1e+05)), silent = TRUE)
      }
    }
    
    # eff_nloptr <- nloptr::nloptr(E0,
    #                              eval_f= f_MP_nloptr,
    #                              lb = efs.min,
    #                              ub = efs.max,
    #                              eval_g_ineq = g_ineq_MP_nloptr ,
    #                              opts = opts,
    #                              q.m = q.m, alpha.m = alpha.m, beta.m = beta.m, pr.m = pr.m,  Cr.f = Cr.f, fc = fc,
    #                              ret.m = ret.m, wd.m = wd.m, wl.m = wl.m, vc.m = vc.m, N = N,  B = B,  K=K,  rho = rho,
    #                              effort.restr = effort.restr, crewS = crewS, catch.restr = catch.restr, tacos = tacos)
    # Et.res[i]   <- sum(eff_nloptr$solution)
    #! if(Et.res[i]>K) Et.res[i] <- K

    
    if(class(eff_opt) != "try-error"){
      res <- K/(1 + exp(-eff_opt[[1]]))
    }else{
      res <- K/(1 + exp(-X))
    }
    Et.res[i] <- sum(res)
    efs.res[, i] <- res/sum(res)

    
    
    # # CHECKS
    # cat('CONVERGENCE: ', eff_opt$convergence, '\n')
    # cat('- capacity cond: ', K >= Et.res, '\n')
    # cat('- efsMax cond  :', names(res) , '\n')
    # cat('                ', res <= efs.max , '\n')
    # cat('- efsMin cond  :', names(res) , '\n')
    # cat('                ', res >= efs.min , '\n')
    # cat('- overshoot    : \n')
    # for (st in stocks.restr) {
    #   Nst <- array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
    #   if(dim(Nst)[1] > 1){
    #     Cam.st <- CobbDouglasAge(sum(res),Nst, wl.m[[st]], wd.m[[st]],
    #                              ret.m[[st]],q.m[[st]],matrix(res/sum(res),ncol = 1),alpha.m[[st]],beta.m[[st]],rho[st])
    #     Cam.st[is.na(Cam.st)] <- 0
    #     if(catch.restr == 'landings') Cam.st <- Cam.st*ret.m[[st]]
    #     Ctot <- ifelse( sum(Cam.st)==0, 1e-08*0.99,sum(Cam.st))
    #   } else{
    #     Cm.st <- CobbDouglasBio(sum(res),Nst, wl.m[[st]], wd.m[[st]],
    #                            q.m[[st]],matrix(res/sum(res),ncol = 1),alpha.m[[st]],beta.m[[st]], ret.m[[st]],rho[st])
    #     if(catch.restr == 'landings') Cm.st <- Cm.st*c(ret.m[[st]])
    #     Ctot <- sum(Cm.st)
    #   }
    #   cat('  * ',st,' -', Ctot <= Cr.f[st], '\n \n')
    # }
    
    if(class(eff_opt) != "try-error"){
      cat("Effort share: ", efs.res[, i], ", ~~~~~ Effort: ", 
        Et.res[i], ", ~~~~~ Funct. Value: ", eff_opt$value, 
        "\n")
    }else{
      cat("Optimization was not successful. Original values used. ", "Effort share: ", efs.res[, i], ", ~~~~~ Effort: ", 
        Et.res[i],   
        "\n")      
    } 
    
    # LO: TB checked!!!!!!
    
    if(LO == TRUE){ # IF LandObl == TRUE => restr == 'catch
      
      lo_res <- MaxProfit_Extra_LO(biols = biols, fleets = fleets, fleets.ctrl = fleets.ctrl, advice.ctrl = advice.ctrl, fl = fl, Et.res = Et.res, 
                                   efs.res = efs.res, efs.min = efs.min, efs.max = efs.max,  yr = yr, ss = ss, 
                                   flnm = flnm, it = it, i = i,  sts = sts, q.m = q.m, alpha.m = alpha.m, 
                                   beta.m = beta.m, pr.m = pr.m,  Cr.f = data.frame(Cr.f), fc = fc,
                                   ret.m = ret.m, wd.m = wd.m, wl.m = wl.m, vc.m = vc.m, N = N,  B = B,  K=K,  rho = rho,
                                   effort.restr = effort.restr, crewS = crewS, catch.restr = restriction, efs.abs = fleets.ctrl[[flnm]]$efs.abs, 
                                   tacos = tacos, tac=TAC[,i, drop=F], Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1,  M = M, Cfyr_1 = Cfyr_1)
    
      
      list2env(lo_res, globalenv())
    } 
  } # END OF ITERATIONS LOOP
  
  
  
  # Update the quota share of this step and the next one if the
  # quota share does not coincide with the actual catch. (update next one only if s < ns).
  
  if(dim(biols[[1]]@n)[4] > 1){ # only for seasonal models
    
    for(st in sts){
      
      
      yr.share       <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [it]
      ss.share       <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, it))# [it,ns]
      quota.share.OR <- matrix(t(yr.share*ss.share), ns, it)
      # The catch.
      catchFun <- fleets.ctrl[[flnm]][[st]][['catch.model']]
      
      Nst  <- N[[st]]
      catchD <- eval(call(catchFun, N = Nst, B = B[st], E = Et.res, efs.m = efs.res, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]], rho = rho[st,]))
      itD <- ifelse(is.null(dim(catchD)), 1, length(dim(catchD)))
      catch <- apply(catchD, itD, sum)  # sum catch along all dimensions except iterations.
      
      quota.share    <- updateQS.SMFB(QS = quota.share.OR, TAC = TAC.yr[st,], catch = catch, season = ss)        # [ns,it]
      
      fleets.ctrl$seasonal.share[[st]][flnm,yr,,,,i] <- t(t(quota.share)/apply(quota.share, 2,sum)) #[ns,it], doble 't' to perform correctly de division between matrix and vector.
    }
  }
  
  
  #  update effort
  fleets[[flnm]]@effort[,yr,,ss] <- Et.res
  for(mt in 1:length(mtnms))  fleets[[flnm]]@metiers[[mt]]@effshare[,yr,,ss] <- efs.res[mt,]
  
  
  return(list(fleets = fleets, fleets.ctrl = fleets.ctrl))
}





#-------------------------------------------------------------------------------
# f_MP_nloptr(E) :: Objective function
#                       (function to be maximized in a fleet by fleet case)
#    - E: numeric(nmt)     sum(E) = Total Effort  across metiers.
#
#    - qa.mt ~ alpha.a.mt ~ beta.a.mt ~ pra.mt.i :: lists with one element per
#                   stock, each element with the form: array(na_st,nmt,it)
#    - Ba ~ list with one element per stock, each element with the form: array(na_st,it)
#
# q.m,...,vc.m must correspond with one iteration of the variable at metier level
# and must have dimension  [nmt,na,nu,1]
# N: alist with one element per stock and each element with dimension [nmt,na,nu]
# Cr.f:[nst,1] quota share
# rho: [nst]
#-------------------------------------------------------------------------------
# 
# f_MP_nloptr <- function(E, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m,
#                         wl.m, N, B, fc, vc.m,   Cr.f,  crewS, K , effort.restr, catch.restr, tacos, rho){
#   
#   nmt <- length(E)
#   
#   res <- 0
#   
#   Cst <- Lst <-  numeric(length(q.m))
#   names(Cst) <- names(Lst) <-names(q.m)
#   
#   #   cat( '**************************************************************************\n')
#   
#   for(st in names(q.m)){
#     
#     #     E1  <- array(E,dim = dim(q.m[[st]]))   # [nmt,na,nu]
#     Nst  <- N[[st]]#array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
#     
#     if(dim(Nst)[1] > 1){
#       Cam <- CobbDouglasAge(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
#                             q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,]) # only 1 iter, MaxProf is applied by iter.
#     }
#     else{
#       Cam <- CobbDouglasBio(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
#                             q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,]) # only 1 iter, MaxProf is applied by iter.
#       Cam <- array(Cam, dim = dim(q.m[[st]]))
#     }
#     
#     Cst[st] <- sum(Cam)
#     Lst[st] <- sum(ret.m[[st]]*Cam)  # multiply the retention vector if landing is the restriction.
#     
#     # The oversized discards are always discarded, but if landing obligation is in place they account in quota (catch == TRUE).
#     if(catch.restr  != 'catch') Cst[st] <- Lst[st]# The restriction is landings.
#     
#     # TAC overshot can be landed or discarded. In the case of landing obligation it is 'discarded' because it does not
#     # contribute to the revenue but it goes against the TAC => TACOS == TRUE
#     if(tacos[st] == FALSE) # TAC Overshot is not discarded.
#       Lrat <- 1
#     else Lrat <- ifelse(Cr.f[st]/Cst[st] > 1, 1, Cr.f[st]/Cst[st])  # TAC Overshot is  discarded.
#     # The overquota discards are proportional to the catch in all the metiers.
#     #   cat('C: ',Cst[st], ', CS: ', Cr.f[st],'\n')
#     res <- res + sum(ret.m[[st]]*Cam*pr.m[[st]])*Lrat
#     #    cat(st, ' - ', sum(ret.m[[st]]*Cam*pr.m[[st]])*Lrat, ' - ', round(Lrat,3), ' - C: ',sum(ret.m[[st]]*Cam)*Lrat,'\n')
#     
#     # cat(st,  ' - ', round(Lrat,3), ' - L: ',sum(ret.m[[st]]*Cam)*Lrat,'\n')
#     
#   }
#   
#   resF <- (1-crewS)*res - sum(vc.m*E) - fc
# 
#   
#   # cat('prof: ', resF, ', rev: ',res, ', creWS:', crewS, ', TCS: ', crewS*res, ', Vc: ', sum(vc.m*E), ', FC: ', fc, '\n' )
#   
#   return(-resF/1e6)
# }
# 


#-------------------------------------------------------------------------------
# f_MP_nloptr(E) :: PENALIZED Objective function
#                       (function to be maximized in a fleet by fleet case)
#    - E: numeric(nmt)     sum(E) = Total Effort  across metiers.
#
#    - qa.mt ~ alpha.a.mt ~ beta.a.mt ~ pra.mt.i :: lists with one element per
#                   stock, each element with the form: array(na_st,nmt,it)
#    - Ba ~ list with one element per stock, each element with the form: array(na_st,it)
#
# q.m,...,vc.m : must correspond with one iteration of the variable at metier level
#                and must have dimension  [nmt,na,nu,1]
# N            : a list with one element per stock and each element with dimension [nmt,na,nu]
# Cr.f         : [nst,1] quota share
# rho          : [nst]
#-------------------------------------------------------------------------------

f_MP_nloptr_penalized <- function(X, efs.min, efs.max, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m,
                                  wl.m, N, B, fc, vc.m,   Cr.f,  crewS, K , effort.restr, catch.restr, stocks.restr, 
                                  efs.abs, tacos, rho, tac, Cyr_1, Nyr_1, Myr_1,  M, Cfyr_1, 
                                  flnm, fleets.ctrl){
  
  E <- K/(1+exp(-X))
  
  nmt <- length(E)
  
  res <- 0
  
  pen_OverShoot <- 0
  
  Cst <- Lst <-  numeric(length(q.m))
  names(Cst) <- names(Lst) <-names(q.m)
  
  resTAC <- numeric(nrow(B))
  names(resTAC) <- rownames(B)
  
  #   cat( '**************************************************************************\n')
  
  for(st in names(q.m)){
    
    #     E1  <- array(E,dim = dim(q.m[[st]]))   # [nmt,na,nu]
    Nst  <- N[[st]] # array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
    
    # if(dim(Nst)[1] > 1){
    #   Cam <- CobbDouglasAge(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
    #                         q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,]) # only 1 iter, MaxProf is applied by iter.
    # }
    # else{
    #   Cam <- CobbDouglasBio(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
    #                         q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,]) # only 1 iter, MaxProf is applied by iter.
    #   Cam <- array(Cam, dim = dim(q.m[[st]]))
    # }
    
    Cam <- eval(call(fleets.ctrl[[flnm]][[st]][['catch.model']], Cr = Cr.f[st,,drop=F], E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
                               q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,],
                               tac=tac[st,,drop=F], Cyr_1 = Cyr_1[[st]], Nyr_1 = Nyr_1[[st]], Myr_1 = Myr_1[[st]],  M = M[[st]], Cfyr_1 = Cfyr_1[[st]],flnm = flnm))
    
    Cam <- array(Cam, dim = dim(q.m[[st]]))
    
    
    Cam[is.na(Cam)] <- 0
    
    Cst[st] <- sum(Cam)
    Lst[st] <- sum(ret.m[[st]]*Cam)  # multiply the retention vector if landing is the restriction.
    
    # The oversized discards are always discarded, but if landing obligation is in place they account in quota (catch == TRUE).
    if(catch.restr  != 'catch') Cst[st] <- Lst[st]# The restriction is landings.
    
    # TAC overshot can be landed or discarded. In the case of landing obligation it is 'discarded' because it does not
    # contribute to the revenue but it goes against the TAC => TACOS == TRUE
    if(tacos[st] == FALSE) # TAC Overshot is not discarded.
      Lrat <- 1
    else Lrat <- ifelse(Cr.f[st,]/Cst[st] > 1, 1, Cr.f[st,]/Cst[st])  # TAC Overshot is  discarded.
    # The overquota discards are proportional to the catch in all the metiers.
    #   cat('C: ',Cst[st], ', CS: ', Cr.f[st],'\n')
    res <- res + sum(ret.m[[st]]*Cam*pr.m[[st]])*Lrat
    #    cat(st, ' - ', sum(ret.m[[st]]*Cam*pr.m[[st]])*Lrat, ' - ', round(Lrat,3), ' - C: ',sum(ret.m[[st]]*Cam)*Lrat,'\n')
    
    # cat(st,  ' - ', round(Lrat,3), ' - L: ',sum(ret.m[[st]]*Cam)*Lrat,'\n')
    
  }
  
  resF <- (1-crewS)*res - sum(vc.m*E) - fc
  
#  cat('income: ', res,', vcost: ', sum(vc.m*E),', crewS: ', crewS*res, ', fcost: ', fc, '\n')
#  cat('profits: ', resF,'effort: ', E,'\n')
  
  #---------------------------------------------------------------------------
  # Incorporate taxes
  #---------------------------------------------------------------------------
  
  taxes <- 0
  
  if (!is.null(fleets.ctrl[[flnm]][[st]][['tax.model']])) {
    
    tax.model <- fleets.ctrl[[flnm]][[st]][['tax.model']]
    
    if (tax.model == "convexTax") {
      
      # Parameters (same units as prices; e.g. eur/ton)
      # - taxes per tonne caught
      gammaC  <- fleets.ctrl[[flnm]][[st]][['gammaC']]
      # - taxes per each tonne that exceeds TAC
      gammaOS <- fleets.ctrl[[flnm]][[st]][['gammaOS']]
      
      for(st in names(q.m))
        taxes <- taxes + 
        gammaC * Cst[st] + gammaOS * ifelse(Cst[st] - Cr.f[st,] < 0, 0, Cst[st] - Cr.f[st,])
      
    } else if (tax.model == "linearTax") {
      
      stop("Linear tax still not available.")
      
    } else if (tax.model == "quadraticTax") {
      
      # Alternative by Helge
      
      # Parameters (same units as prices; e.g. eur/ton)
      # - taxes per tonne caught
      gammaC  <- fleets.ctrl[[flnm]][[st]][['gammaC']]
      # - taxes per each tonne that exceeds TAC
      gammaOS <- fleets.ctrl[[flnm]][[st]][['gammaOS']]
      
      for(st in names(q.m))
        taxes <- taxes + 
          gammaC * Cst[st]/Cr.f[st,] + gammaOS/2 * Cr.f[st,] * (Cst[st]/Cr.f[st,])^2
      
    }
    
  } else
    
  
  resF <- resF - taxes
  
  
  #---------------------------------------------------------------------------
  # constraint on effort-share: absolute or relative values.
  #---------------------------------------------------------------------------
  
  if(efs.abs == TRUE){
    pen_efsMax <- sum(log(E/(efs.max - E)))
    # Using the inverse 1/efs and using efs.min, we bound the minimum
    pen_efsMin <- sum(log((1/E)/((1/efs.min) - (1/E))))
  }
  else{
    efs <- E/sum(E)
    pen_efsMax <- pen_efsMin <- 0
    #efs <- efs#to allow 0 values in efs
    if(length(efs) > 1){ # If there is onlly one metier we don't want any penalty on this.
      # This penalty ensures that efs < efs.max because otherwise  log(efs/(efs.max - efs)) = NaN
      pen_efsMax <- sum(log(efs/(efs.max - efs)))
      # Using the inverse 1/efs and using efs.min, we bound the minimum
      pen_efsMin <- sum(log((1/efs)/((1/efs.min) - (1/efs))))
    }}
  
  
  # print(efs)
  # print(efs.min)
  
  #---------------------------------------------------------------------------
  # constraint on capacity
  #---------------------------------------------------------------------------
  Et <- sum(E)
  Et <- Et  #to allow 0 values in efs
  pen_Et <- log(Et/(K-Et))
  
  
  #---------------------------------------------------------------------------
  # constraint on catches, comply with all the TACS ('min') or only with one.
  #---------------------------------------------------------------------------
  if(effort.restr == 'none') resTAC <- 0 # The overshoot of TACs is not penalized but it does not produce any rent so overshooting 
  # all of them should not be wroth unless there are some unconstrained stocks
  # like OTH, in this case if no TAC and the biomas sin unlimited, there must be some stock constraint.  
  if(effort.restr == 'min'){
    
    resTAC <- rep(0, length(q.m))  # One resctriction per stock.
    names(resTAC) <- names(q.m)
    
    for(st in names(q.m)){
      Nst <- N[[st]] #array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
      if(dim(Nst)[1] > 1){
        Cam <- CobbDouglasAge(sum(E),Nst, wl.m[[st]], wd.m[[st]],
                              ret.m[[st]],q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]],rho[st,])
        
        Cam[is.na(Cam)] <- 0
        
        if(catch.restr == 'landings') Cam <- Cam*ret.m[[st]]
        
        #       resTAC[st] <- 1/(1+2^((-Cr.f[st]+sum(Cam))/0.00005)) this constraint does not work for all the stocks simultaneously
        Ctot <- ifelse( sum(Cam)==0, 1e-08*0.99,sum(Cam))
        resTAC[st] <- log(Ctot/(Cr.f[st,]-Ctot))
        #   cat(st, ' - ', sum(Cam), '\n')
      }
      else{
        Cm <- CobbDouglasBio(sum(E),Nst, wl.m[[st]], wd.m[[st]],
                             q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]], ret.m[[st]],rho[st,])
        
        if(catch.restr == 'landings') Cm <- Cm*c(ret.m[[st]])
        
        # resTAC[st] <- sum(Cm)  - Cr.f[st]
        #       cat(st, ' - ', sum(Cm), '\n')
        #
        # resTAC[st] <- 1/(1+2^((-Cr.f[st]+sum(Cm))/0.00005)) this constraint does not work for all the stocks simultaneously
        Ctot <- ifelse( sum(Cm)==0, 1e-08*0.99,sum(Cm))
        resTAC[st] <- log(Ctot/(Cr.f[st,]-Ctot))
        
      }
      
    }
    pen_OverShoot <- sum(resTAC[stocks.restr]) # only the overshoot of the TAC of some stocks penalizes the function.
  }
  
  if(!(effort.restr %in% c('none', 'min'))){
    stk.cnst <- effort.restr
    
    Nst <- N[[stk.cnst]] # array(N[[stk.cnst]][drop=T],dim = dim(N[[stk.cnst]])[c(1,3,6)])
    

    Cm <- eval(call(fleets.ctrl[[flnm]][[stk.cnst]][['catch.model']], Cr = Cr.f[stk.cnst,,drop=F], E = sum(E), N = Nst, wl.m = wl.m[[stk.cnst]], wd.m = wd.m[[stk.cnst]], ret.m = ret.m[[stk.cnst]],
                       q.m = q.m[[stk.cnst]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[stk.cnst]], beta.m = beta.m[[stk.cnst]], rho = rho[stk.cnst,],
                       tac=tac[stk.cnst,], Cyr_1 = Cyr_1[[stk.cnst]], Nyr_1 = Nyr_1[[stk.cnst]], Myr_1 = Myr_1[[stk.cnst]],  M = M[[stk.cnst]], Cfyr_1 = Cfyr_1[[stk.cnst]],flnm = flnm))
 
    Cm <- array(Cm, dim = dim(wl.m[[stk.cnst]]))
    
    
    if(catch.restr == 'landings') Cm <- Cm*c(ret.m[[stk.cnst]])
    
    # resTAC[st] <- 1e7*sum(1/(1+2^((-Cr.f[stk.cnst]+sum(Cm))/0.00005))) This penalty works properly only if the multiplier is in the scale of the profits.
    resTAC[stk.cnst] <- log(sum(Cm)/(Cr.f[stk.cnst,]-sum(Cm)))
    
    pen_OverShoot <- resTAC[stk.cnst]
  }
  
  # cat('-----------------------------------------\n')
  # cat('profits: ', resF, '\n')
  # cat('overshoot: ', pen_OverShoot, '\n')
  # cat('efsMax: ', pen_efsMax, '\n')
  # cat('efsMin: ', pen_efsMin, '\n')
  
  # cat('prof: ', resF, ', rev: ',res, ', creWS:', crewS, ', TCS: ', crewS*res, ', Vc: ', sum(vc.m*E), ', FC: ', fc, '\n' )
  
  obj <-  -sum(resF) - pen_OverShoot + pen_efsMax + pen_efsMin + pen_Et
  
  # cat('***** Obj ********: ', obj, '\n')
  
  return(obj/1e6)
}


#-------------------------------------------------------------------------------
# nlin.maxprofits(E)
#   * The maximization of profits is restringed to the compliance,
#     in some degree, of the TACs. There will be different options based on the
#     Fcube like approach. The contraint for all the stock will have similar form
#     so we write the general function  'nlin.maxprofits' and then we will use
#     it in accordance with the option choosen.
# MIN OPTION
#   * The maximization of benefits is constrained to the TACs of all the stocks
#     (if one stock is not managed put each TAC share equal to infinity).
#     That is, the catch must be below TAC quota share.
#-------------------------------------------------------------------------------
# g_ineq_MP_nloptr <- function(E, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m,
#                              wl.m, N, B, fc, vc.m,  Cr.f,  crewS, K, effort.restr, catch.restr, tacos, rho){
#   
#   nmt <- length(E)
#   
#   stnms <- names(N)
#   
#   res <- 0
#   
#   Cst <- Lst <-  NULL
#   #   cat( '**************************************************************************\n')
#   # constraint on catches, comply with all the TACS ('min') or only with one.
#   if(effort.restr == 'min'){
#     
#     resTAC <- rep(0, length(q.m))  # One resctriction per stock.
#     names(resTAC) <- names(q.m)
#     
#     for(st in names(q.m)){
#       Nst <- N[[st]] #array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
#       if(dim(Nst)[1] > 1){
#         Cam <- CobbDouglasAge(sum(E),Nst, wl.m[[st]], wd.m[[st]],
#                               ret.m[[st]],q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]],rho[st,])
#         
#         if(catch.restr == 'landings') Cam <- Cam*ret.m[[st]]
#         
#         resTAC[st] <- sum(Cam)  - Cr.f[st]
#         #   cat(st, ' - ', sum(Cam), '\n')
#       }
#       else{
#         Cm <- CobbDouglasBio(sum(E),Nst, wl.m[[st]], wd.m[[st]],
#                              q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]], ret.m[[st]],rho[st,])
#         
#         if(catch.restr == 'landings') Cm <- Cm*c(ret.m[[st]])
#         
#         resTAC[st] <- sum(Cm)  - Cr.f[st,]
#         #       cat(st, ' - ', sum(Cm), '\n')
#       }
#     }}
#   else{
#     stk.cnst <- effort.restr
#     
#     Nst <- N[[stk.cnst]]        #array(N[[stk.cnst]][drop=T],dim = dim(N[[stk.cnst]])[c(1,3,6)])
#     
#     if(dim(Nst)[1] > 1){
#       Cm <- CobbDouglasAge(sum(E),Nst, wl.m[[stk.cnst]], wd.m[[stk.cnst]],
#                            ret.m[[stk.cnst]],q.m[[stk.cnst]],matrix(E/sum(E),ncol = 1),alpha.m[[stk.cnst]],beta.m[[stk.cnst]],rho[stk.cnst,])   
#     }
#     else{
#       Cm <- array(CobbDouglasBio(sum(E),Nst, wl.m[[stk.cnst]], wd.m[[stk.cnst]], q.m[[stk.cnst]],matrix(E/sum(E),ncol = 1),alpha.m[[stk.cnst]],beta.m[[stk.cnst]], ret.m[[stk.cnst]],rho[stk.cnst,]),
#                   dim = dim(wl.m[[stk.cnst]]))
#     }
#     
#     if(catch.restr == 'landings') Cm <- Cm*c(ret.m[[stk.cnst]])
#     
#     resTAC <- sum(Cm) - Cr.f[stk.cnst,]
#   }
#   
#   # constraint on capacity.
#   resK <- sum(E)-K
#   
#   return(c(resTAC, resK))
# }
# 

#******************************************************************************************
#-------------------------------------------------------------------------------
#  FUNCTION FOR OPTIMIZATION IN LANDING OBLIGATION SCENARIOS
#
#-------------------------------------------------------------------------------
#******************************************************************************************


#-------------------------------------------------------------------------------
# f_MP_LO_nloptr(X) :: Objective function, X = v(E,tau)
#                       (function to be maximized in a fleet by fleet case)
#    - E: numeric(nmt)     sum(E) = Total Effort  across metiers.
#    - tau: the percentage of quota swaped to restrictive stock.
#    - qa.mt ~ alpha.a.mt ~ beta.a.mt ~ pra.mt.i :: lists with one element per
#                   stock, each element with the form: array(na_st,nmt,it)
#    - Ba ~ list with one element per stock, each element with the form: array(na_st,it)
#
# q.m,...,vc.m must correspond with one iteration of the variable at metier level
# and must have dimension  [nmt,na,nu,1]
# N: alist with one element per stock and each element with dimension [nmt,na,nu]
# Cr.f:[nst,1] quota share
# rho: [nst]
# 
#  Effort restriction is always 'min', i.e all the Quotas must be fulfiled.
#-------------------------------------------------------------------------------

f_MP_LO_nloptr <- function(X, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m,
                           wl.m, N, B, fc, vc.m,   Cr.f,  crewS, K, rho,STRs,STDs, Cr.f.new, tau.old, lambda.lim){


  E <- X[-(1:length(STRs))]

  tau     <- X[1:length(STRs)]

  nmt <- length(E)

  res <- 0

  Cst <- Lst <-  numeric(length(q.m))
  names(Cst) <- names(Lst) <-names(q.m)

  #   cat( '**************************************************************************\n')

  for(st in names(q.m)){

    #     E1  <- array(E,dim = dim(q.m[[st]]))   # [nmt,na,nu]
    Nst  <- N[[st]]  # array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])

    if(dim(Nst)[1] > 1){
      Cam <- CobbDouglasAge(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
                            q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,]) # only 1 iter, MaxProf is applied by iter.
    }
    else{
      Cam <- CobbDouglasBio(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
                            q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st,]) # only 1 iter, MaxProf is applied by iter.
      Cam <- array(Cam, dim = dim(q.m[[st]]))
    }

    Cst[st] <- sum(Cam) # IN LO 'CATCH' IS **ALWAYS** THE RESTRICTOR
    #    Lst[st] <- sum(ret.m[[st]]*Cam)  # multiply the retention vector if landing is the restriction.

    # The oversized discards are always discarded, but if landing obligation is in place they account in quota (catch == TRUE).
    #   if(catch.restr  != 'catch') Cst[st] <- Lst[st]# The restriction is landings.

    # TAC overshot can be landed or discarded. In the case of landing obligation it is 'discarded' because it does not
    # contribute to the revenue but it goes against the TAC => TACOS == TRUE
    # IN LO **EVERYTHING** IS LANDED.
    #   if(tacos[st] == FALSE) # TAC Overshot is not discarded.
    Lrat <- 1
    #   else Lrat <- ifelse(Cr.f[st]/Cst[st] > 1, 1, Cr.f[st]/Cst[st])  # TAC Overshot is  discarded.
    # The overquota discards are proportional to the catch in all the metiers.
    res <- res + sum(ret.m[[st]]*Cam*pr.m[[st]])*Lrat


  }

  resF <- (1-crewS)*res - sum(vc.m*E) - fc

  # cat('prof: ', resF, ', rev: ',res, ', creWS:', crewS, ', TCS: ', crewS*res, ', Vc: ', sum(vc.m*E), ', FC: ', fc, '\n' )

  return(-resF/1e6)
}

#-------------------------------------------------------------------------------
# Constraints in the swap of quotas in MP.
#   o Total effort in each metier is bigger than the a certain level. (previous effort in the algorithm)
#   o The catch of the stocks != 'str' and 'std' is lower than the 'new catch quota (catch in Cr.f.new).
#   o tau_old < tau_new < 0.9 and:
#       - Cstr < Qstr + (tau_old + tau_new)*Qstd
#       - Cstd < QNEWstd - tau_new*Qstd  
#  Effort restriction is always 'min', i.e all the Quotas must be fulfiled.
#-------------------------------------------------------------------------------
g_ineq_MP_LO_nloptr <- function(X, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m, wl.m, N, B, fc, vc.m,
                                Cr.f, crewS, K, rho, STRs,STDs, Cr.f.new, tau.old, lambda.lim){
  E <- X[-(1:length(STRs))]

  tau <- X[1:length(STRs)]
  names(tau) <- STRs

  # browser()

  nmt <- length(E)

  stnms <- names(N)

  res <- 0

  Cst <- Lst <-  NULL
  #   cat( '**************************************************************************\n')
  # constraint on catches, comply with all the TACS ('min') or only with one.

  resTAC <- rep(0, length(q.m))  # One resctriction per stock.
  names(resTAC) <- names(q.m)

  # NEW quotas for str and std
  for(str in STRs){
    Cr.f.new[str] <- Cr.f[str] + tau[str]*Cr.f[STDs[,str]]
    Cr.f.new[STDs[,str]] <- Cr.f.new[STDs[,str]] - (tau[str] - tau.old[str])*Cr.f[STDs[,str]]
  }

  for(st in names(q.m)){
    Nst <- N[[st]] #array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
    if(dim(Nst)[1] > 1){
      Cam <- CobbDouglasAge(sum(E),Nst, wl.m[[st]], wd.m[[st]],
                            ret.m[[st]],q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]],rho[st,])

      resTAC[st] <- sum(Cam)  - Cr.f.new[st]
      #   cat(st, ' - ', sum(Cam), '\n')
    }
    else{
      Cm <- CobbDouglasBio(sum(E),Nst, wl.m[[st]], wd.m[[st]],
                           q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]], ret.m[[st]],rho[st,])
      resTAC[st] <- sum(Cm)  - Cr.f.new[st]
      #       cat(st, ' - ', sum(Cm), '\n')
    }

  }
  #  browser()
  lambda <- -lambda.lim

  for(st in unique(STDs[1,])){
    strs_std       <- names(STDs[1,which(STDs[1,] == st)])
    lambda[st] <- lambda[st] + sum(tau[strs_std] - tau.old[strs_std])
  }

  resK <- sum(E)-K
  # print(c(resTAC, resK, lambda))
  # return(c(lambda))
  return(c(resTAC, resK, lambda))
}

