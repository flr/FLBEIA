#-------------------------------------------------------------------------------
# Functions to be used within Rdonlp2 to obtain effort and effortshare according
# to maximization of profits, restringed to Fcube like approach in 
# TAC constraints.
#
# Dorleta Garcia - Azti Tecnalia
# Created: 22/11/2011 19:02:56
# Changed: 23/11/2011 10:16:01
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# MaxProfit.stkCnst(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1)
# Catch production function based on Cobb-Doug at age level.
# OVER-QUOTA LANDINGS ARE DISCARDED => NOT BENEFIT FROM THEM.
#-------------------------------------------------------------------------------
MaxProfit <- function(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1){
    
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
    
    # Check fleets.ctrl elements.
 # This checks the restriction for all the fleets and we only need to check the fleet we are proyecting.
 #   if(! all(sapply(names(fleets), function(x) fleets.ctrl[[x]]$restriction %in% c('catch', 'landings'))))
 #       stop("fleets.ctrl$restriction must be equal to 'catch' or 'landings'")
    if(!(fleets.ctrl[[flnm]]$restriction %in% c('catch', 'landings') )) stop("fleets.ctrl$restriction for fleet, ', flnm, ', must be equal to 'catch' or 'landings'")

    # Dimensions.
    nst <- length(biols);          stnms <- names(biols)
    ns  <- dim(biols[[1]]@n)[4]
    it  <- dim(biols[[1]]@n)[6]
    flnms <- names(fleets)
    
    nmt <- length(fleets[[flnm]]@metiers)

    Et.res  <- numeric(it)
    efs.res <- matrix(NA,nmt,it)


    for(i in 1:it){
        # Biomass at age.
        B    <- sapply(stnms, function(x){   # biomass in the middle of the season  [nst]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss,,i, drop=T])
                                else return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss,,i, drop=T])})

        N   <- lapply(stnms, function(x){   # biomass at age in the middle  of the season, list elements: [na,1,nu,1,1,1]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss,,i, drop = FALSE])
                                else return((biols[[x]]@n)[,yr,,ss,,i, drop = F])})
        names(N) <- stnms

        QS.fls   <- sapply(stnms, function(x){           # matrix [nf,nst]
                            # Calculate QS by fleet for the year and season
                            yr.share    <- advice$quota.share[[x]][,yr,,,,i, drop=T]        # [nf]
                            ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss,,i, drop=T]   # [nf]
                            quota.share <- yr.share*ss.share # [nf]
                            quota.share[is.na(quota.share)] <- 0
                            return(quota.share)})
         names(QS.fls) <- stnms

         # If TAC >= B*alpha => TAC = B*alpha.
         TAC.yr   <- advice$TAC[,yr,,,,i,drop=T]    # nst
         rho       <- fleets.ctrl$catch.threshold[,yr,,ss,,i, drop=T]  # [ns]
         QS.ss    <- colSums(QS.fls)  # [nst] Total seasonal quota share

         TAC <- ifelse(B*rho < TAC.yr*QS.ss, B*rho, TAC.yr*QS.ss)

         # Re-scale QS to fleet share within the season instead of season-fleet share within year.
         QS   <- sweep(QS.fls, 2, apply(QS.fls,2, sum), "/")  # [nf,nst]
         QS[is.na(QS)] <- 0

         fl    <- fleets[[flnm]]
         sts   <- catchNames(fl)
         mtnms <- names(fl@metiers)
         nmt   <- length(mtnms)

         # flinfo: matrix with information on which metier catch which stock.
         fl.        <- FLFleetsExt(fl)
         names(fl.) <- flnm
         flinfo     <- stock.fleetInfo(fl.)
         flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')

         efs.m <- sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr-1,,ss,,i, drop=T]) #[nmt], previous year effort share because we don't know the effort share this year.

         vc.m <- sapply(mtnms, function(x) fl@metiers[[x]]@vcost[,yr,,ss,,i, drop=T])  #[nmt]

         fc    <- fl@fcost[,yr,,ss, drop=T]*covars$NumbVessels[flnm,yr,,ss,,i, drop=T] # [1]

         crewS <- fl@crewshare[,yr,,ss,,i, drop=T] # [i]

         effs <- numeric(nst); names(effs) <- stnms
         Cr.f <- numeric(nst); names(Cr.f) <- stnms

         q.m <- alpha.m <- beta.m  <- pr.m <- ret.m <- wd.m <- wl.m <-vector('list', length(sts))
         names(q.m) <- names(pr.m) <- names(alpha.m) <- names(beta.m) <- names(ret.m) <- names(wl.m) <- names(wd.m) <- sts

         tacos <- logical(length(sts))
         names(tacos) <- sts

         for(st in sts){     # q.m, alpha.m.... by metier but stock specific

                # identify the first metier that catch stock st
                mtst <- flinfo[[st]][2]

                age.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[1]]
                age.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[1]]
                age.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[1]]
                age.pr    <- dimnames(fl@metiers[[mtst]]@catches[[st]]@price)[[1]]

                unit.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[3]]
                unit.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[3]]
                unit.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[3]]
                unit.pr    <- dimnames(fl@metiers[[mtst]]@catches[[st]]@price)[[3]]

                q.m[[st]]     <- array(0, dim = c(length(mtnms), length(age.q),     length(unit.q), 1),      dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1))
                alpha.m[[st]] <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), 1), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1))
                beta.m[[st]]  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
                ret.m[[st]]   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
                wl.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
                wd.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
                pr.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))


                # if TAC overshoot is not discarded, it is sold and it contributes to the revenue.
                tacos[st] <- ifelse(is.null(fleets.ctrl[[flnm]][[st]][['discard.TAC.OS']]), TRUE,fleets.ctrl[[flnm]][[st]][['discard.TAC.OS']])

                for(mt in mtnms){
                       if(!(st %in% names(fl@metiers[[mt]]@catches))) next

                       q.m[[st]][mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, ,i,drop = TRUE]
                       alpha.m[[st]][mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, ,i,drop = TRUE]
                       beta.m[[st]][mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, ,i,drop = TRUE]
                       ret.m[[st]][mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, ,i,drop = TRUE]
                       wl.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, ,i,drop = TRUE]
                       wd.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, ,i,drop = TRUE]

                       # The price is taken from the year before, because price for the current year is updated after catch is produced,
                       # if the price was dynamically updated inside this function the optimizer could crash.
                       pr.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@price[,yr-1,,ss,,i, drop = TRUE]
                       }
                       Cr.f[st] <- TAC[st]*QS[flnm,st]
         }

        
         # Calculate the starting point based on the effort that correspond with the TAC quotas.
        effs <- numeric(length(q.m))
        names(effs) <- names(q.m)
        
        for(st in names(q.m)){
            effort.fun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'effort', sep = '.')
            Nst  <- array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
            effs[st] <-  eval(call(effort.fun, Cr = Cr.f[st],  N = Nst[,,i,drop=F], q.m = q.m[[st]][,,,i,drop=F],
                                  efs.m = matrix(efs.m,nmt,1), alpha.m = alpha.m[[st]][,,,i,drop=F], beta.m = beta.m[[st]][,,,i,drop=F],
                                  ret.m = ret.m[[st]][,,,i,drop=F], wl.m = wl.m[[st]][,,,i,drop=F], wd.m = wd.m[[st]][,,,i,drop=F],
                                  restriction = fleets.ctrl[[flnm]]$restriction))
        }
        
        effort.restr <- fleets.ctrl[[flnm]]$effort.restr
        K <- c(fl@capacity[,yr,,ss,,it,drop=T])
        
        Et  <- 0.9*ifelse(effort.restr == 'min', min(effs), effs[effort.restr]) 
        Et <- ifelse(Et < K, Et, K*0.9)

  #       browser()

         catch.rest <- ifelse(is.null(fleets.ctrl[[flnm]]$restriction), 'landings', fleets.ctrl[[flnm]]$restriction)

         eff_nloptr <- nloptr(Et*efs.m,
             eval_f= f_MP_nloptr,
             lb = rep(0, nmt),
             ub = rep(K, nmt),
             eval_g_ineq = g_ineq_MP_nloptr ,
             opts = list("algorithm" = "NLOPT_LN_COBYLA", maxeval = 1e9, xtol_abs = rep(1e-4,nmt), xtol_rel = 1e-4, maxtime = 300),
             q.m = q.m, alpha.m = alpha.m, beta.m = beta.m, pr.m = pr.m,  Cr.f = Cr.f, fc = fc,
             ret.m = ret.m, wd.m = wd.m, wl.m = wl.m, vc.m = vc.m, N = N,  B = B,  K=K,  rho = rho,
             effort.restr = effort.restr, crewS = crewS, catch.rest = catch.rest, tacos = tacos)
             Et.res[i]   <- sum(eff_nloptr$solution)
       efs.res[,i] <- eff_nloptr$solution/sum(eff_nloptr$solution)
        cat('Effort share: ', efs.res[,i], ', ~~~~~ Effort: ',Et.res[i], ', ~~~~~ Benefit: ', eff_nloptr$objective, '\n')

       # Update the quota share of this step and the next one if the
       # quota share does not coincide with the actual catch. (update next one only if s < ns).
       for(st in sts){
       #   browser()
       #  if(st == 'SKH') browser()
        
            yr.share       <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [it]
            ss.share       <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, it))# [it,ns]
            quota.share.OR <- matrix(t(yr.share*ss.share), ns, it)
            # The catch.
            catchFun <- fleets.ctrl[[flnm]][[st]][['catch.model']]

            Nst  <- array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
            catchD <- eval(call(catchFun, N = Nst, B = B[st], E = Et.res, efs.m = efs.res, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]], rho = rho[st]))
            itD <- ifelse(is.null(dim(catchD)), 1, length(dim(catchD)))
            catch <- apply(catchD, itD, sum)  # sum catch along all dimensions except iterations.

            quota.share    <- updateQS.SMFB(QS = quota.share.OR, TAC = TAC.yr[st], catch = catch, season = ss)        # [ns,it]
                              
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
#    - E: numeric(nmt)     E[-1] = Tiotal Effort  across metiers.
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

f_MP_nloptr <- function(E, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m,
                wl.m, N, B, fc, vc.m,   Cr.f,  crewS, K , effort.restr, catch.rest, tacos, rho){

    nmt <- length(E)

    res <- 0

    Cst <- Lst <-  numeric(length(q.m))
    names(Cst) <- names(Lst) <-names(q.m)

    for(st in names(q.m)){

   #     E1  <- array(E,dim = dim(q.m[[st]]))   # [nmt,na,nu]
        Nst  <- array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])

        if(dim(Nst)[1] > 1){
            Cam <- CobbDouglasAge(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
                   q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st]) # only 1 iter, MaxProf is applied by iter.
        }
        else{
            Cam <- CobbDouglasBio(E = sum(E), N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]],
                   q.m = q.m[[st]], efs.m = matrix(E/sum(E),ncol = 1), alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st]) # only 1 iter, MaxProf is applied by iter.
            Cam <- array(Cam, dim = dim(q.m[[st]]))
         }

        Cst[st] <- sum(Cam)
        Lst[st] <- sum(ret.m[[st]]*Cam)  # multiply the retention vector if landing is the restriction.

        # The oversized discards are always discarded, but if landing obligation is in place they account in quota (catch == TRUE).
        if(catch.rest  != 'catch') Cst[st] <- Lst[st]# The restriction is landings.

        # TAC overshot can be landed or discarded. In the case of landing obligation it is 'discarded' because it does not
        # contribute to the revenue but it goes against the TAC => TACOS == TRUE
        if(tacos[st] == FALSE) # TAC Overshot is not discarded.
            Lrat <- 1
        else Lrat <- ifelse(Cr.f[st]/Cst[st] > 1, 1, Cr.f[st]/Cst[st])  # TAC Overshot is  discarded.
                                                                              # The overquota discards are proportional to the catch in all the metiers.

        res <- res + sum(ret.m[[st]]*Cam*pr.m[[st]])*Lrat


    }

    res <- (1-crewS)*res - sum(vc.m*E) - fc

    return(-res/1e6)
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
g_ineq_MP_nloptr <- function(E, q.m, alpha.m, beta.m, pr.m, ret.m, wd.m,
                wl.m, N, B, fc, vc.m,  Cr.f,  crewS, K, effort.restr, catch.rest, tacos, rho){

    nmt <- length(E)

    stnms <- names(N)

    res <- 0

    Cst <- Lst <-  NULL

    # constraint on catches, comply with all the TACS ('min') or only with one.
    if(effort.restr == 'min'){

        resTAC <- rep(0, length(q.m))  # One resctriction per stock.
        names(resTAC) <- names(q.m)

        for(st in names(q.m)){
            Nst <- array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
            if(dim(Nst)[1] > 1){
                resTAC[st] <- sum(CobbDouglasAge(sum(E),Nst, wl.m[[st]], wd.m[[st]],
                       ret.m[[st]],q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]],rho[st]))  - Cr.f[st]
            }
            else{
                resTAC[st] <- sum(CobbDouglasBio(sum(E),Nst, wl.m[[st]], wd.m[[st]],
                      q.m[[st]],matrix(E/sum(E),ncol = 1),alpha.m[[st]],beta.m[[st]], ret.m[[st]],rho[st]))  - Cr.f[st]
           }
    }}
    else{
         stk.cnst <- effort.restr

         Nst <- array(N[[stk.cnst]][drop=T],dim = dim(N[[stk.cnst]])[c(1,3,6)])
         resTAC <- sum(CobbDouglasAge(sum(E),Nst, wl.m[[stk.cnst]], wd.m[[stk.cnst]],
                       ret.m[[stk.cnst]],q.m[[stk.cnst]],matrix(E/sum(E),ncol = 1),alpha.m[[stk.cnst]],beta.m[[stk.cnst]],rho[stk.cnst])) - Cr.f[stk.cnst]}

    # constraint on capacity.
    resK <- sum(E)-K

    # constraint on catch, Cst < alpha_st*Bst for any st
#    resB <- numeric(length(stnms))    This constraint makes no sense.
#    names(resB) <- stnms
#    for(st in names(q.m)){
#
#            E1  <- array(E, dim(q.m[[st]]))
#            Nst  <- array(rep(N[[st]][drop=T], each = nmt),dim = dim(q.m[[st]]))
#            resB[st] <- sum(q.m[[st]]*(E1^alpha.m[[st]])*(Nst*(ret.m[[st]]*wl.m[[st]] +
#                    (1-ret.m[[st]])*wd.m[[st]]))^beta.m[[st]])  - B[st]*rho[st]
#    }
 #   resB[] <- -1
  #  resTAC[] <- -1
  #  resK <- -1

    return(c(resTAC, resK))
}


