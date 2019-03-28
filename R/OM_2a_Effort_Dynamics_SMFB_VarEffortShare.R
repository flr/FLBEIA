#-------------------------------------------------------------------------------
#        ORIGINAL SMFB function with extra functionalities to reproduce the 
#       landing obligation policy with the following exemptions:
#         o Minimise
#         o Quota transfers between years.
#         o Quota swap btw stocks.
#
# New arguments in fleets.ctrl[[flnm]] object to control the landing obligation implementation:
#   o LandObl: Logical TRUE/FALSE, Is the landing obligation in place?
#   o LandObl_minimis: logical[nyr], is minimis exemption applied? one element per year with ALL the years, including historical ones.
#   o LandObl_yearTransfer: logical[nyr], is quota transfer between years exemption applied? one element per year with ALL the years, including historical ones.
#   o LandObl_minimis_p: matrix[st,ny], if minimis applied declare the proportion for each year. 
#   o LandObl_yearTransfer_p: matrix[st,ny], if minimis applied declare the proportion for each year.  
#   o LandObl_discount_yrtransfer: If yearTransfer == TRUE, the discount to be applied in the first year. 
#   o LO_stk_grp: The groups to swap quotas. 
#
# 'SMFB' - (Simple mixed fisheries behaviour). - Everything constant  except effort 
#       that is updated based on landings or catch share. 
#       (multiple TACs so min, max Effort options are applied)
#
# Dorleta GarcYYYa
# Created: 23/10/2014 12:33:04
# Changed: 13/01/2015
# Changed: 01/04/2015 Itsaso Carmona 
# Changed: 29/04/2015 Itsaso carmona (LO in some years)
# Added Effort share models: 20/03/2019 Dorleta 
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SMFB_LO(fleets, biols, covars, fleets.ctrl, year = 1, season = 1)
#-------------------------------------------------------------------------------
SMFB_ES <- function(fleets, biols, BDs, covars, advice, biols.ctrl, fleets.ctrl, advice.ctrl, flnm, year = 1, season = 1,...){
    
    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )

    # If year/season/iter numerics => indicate position 
    # else names => get positions.
    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    dimnms <- dimnames(biols[[1]]@n) 
    
    # 'year' dimension.
    yr <- year
    if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
    if(length(yr) == 0) stop('The year is outside object time range')  
    
    # 'season' dimension.
    ss <- season
    if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
    if(length(ss) == 0) stop('The season is outside object season range')  
    
    # Check fleets.ctrl elements.
    restriction <- ifelse(length(fleets.ctrl[[flnm]]$restriction) == 1, fleets.ctrl[[flnm]]$restriction, fleets.ctrl[[flnm]]$restriction[year])
    if(!(restriction %in% c('catch', 'landings')))
        stop("fleets.ctrl[[f]]$restriction must be equal to 'catch' or 'landings'")
     

    # Dimensions.
    stnms <- catchNames(fleets[[flnm]])
    nst <- length(stnms)          
    ns  <- dim(biols[[1]]@n)[4]
    it  <- dim(biols[[1]]@n)[6]
    flnms <- names(fleets)
    
    # Data
    B    <- matrix(t(sapply(stnms, function(x){   # biomass in the middle of the season  [nst,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss, drop=T])
                                else{
                                  if(biols.ctrl[[x]][['growth.model']] == 'fixedPopulation'){
                                    return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss, drop=T])
                                  }
                                  else{
                                    return((biols[[x]]@n*biols[[x]]@wt + BDs[[x]]@gB)[,yr,,ss, drop=T])
                                  }
                                  
                                } })) , nst,it, dimnames = list(stnms, 1:it))

    N   <- lapply(stnms, function(x){   # biomass at age in the middle  of the season, list elements: [na,1,nu,1,1,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss, drop = FALSE])
                                else{
                                  if(biols.ctrl[[x]] == 'fixedPopulation'){
                                    return((biols[[x]]@n)[,yr,,ss, drop=F])
                                  }
                                  else{
                                    return((biols[[x]]@n + BDs[[x]]@gB)[,yr,,ss, drop=F])
                                  } } })
    
    names(N) <- stnms
    

                      
    # Quota share          
    QS   <- lapply(stnms, function(x){           # list of stocks, each stock [nf,it]
                            # Calculate QS by fleet for the year and season
                            yr.share    <- advice$quota.share[[x]][,yr,, drop=T]        # [nf,it]
                            ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss, drop=T]   # [nf,it]
                            quota.share <-  matrix(yr.share*ss.share, length(flnms), it, dimnames = list(flnms, 1:it))
                            quota.share[is.na(quota.share)] <- 0
                            return(quota.share)})         
    names(QS) <- stnms
        
    # If TAC >= B*alpha => TAC = B*alpha.
    TAC.yr  <- matrix(advice$TAC[stnms,yr,drop=T], nst, it, dimnames = list(stnms, 1:it))   # [nst,it]
    rho       <- fleets.ctrl$catch.threshold[,yr,,ss, drop=T]  # [ns,it]
    
    # if rho is a numeric => there is only one stock and one iteration wihout names => name it
    if(is.null(dim(rho)) & is.null(names(rho))) names(rho) <- stnms
      
    QS.ss    <- matrix(t(sapply(stnms, function(x) apply(QS[[x]],2,sum))), nst,it, dimnames = list(stnms, 1:it))  # [nst,it]
                            
    for(stknm in  stnms){
        tacos.fun <- fleets.ctrl[[flnm]][[stknm]]$TAC.OS.model
        if(is.null(tacos.fun))   alpha <- rep(1,it)
        else{
            alpha <- eval(call(tacos.fun, fleets = fleets, TAC = TAC.yr, fleets.ctrl = fleets.ctrl, flnm = flnm, stknm = stknm, year = year, season = season))
        }
        TAC.yr[stknm,] <- TAC.yr[stknm,]*alpha 

    }
    
    if(it > 1){    
      if(length(stnms) == 1) rho <- matrix(rho, 1,it, dimnames = list(stnms, 1:it))
      
      TAC <- ifelse(B*rho[stnms,] < TAC.yr*QS.ss, B*rho[stnms,], TAC.yr*QS.ss)
    }
    else TAC <- ifelse(B*rho[stnms] < TAC.yr*QS.ss, B*rho[stnms], TAC.yr*QS.ss)

    # Re-scale QS to fleet share within the season instead of season-fleet share within year.
    QS   <- lapply(stnms, function(x){          # list of stocks, each stock [nf,it]
                            res <- sweep(QS[[x]], 2, apply(QS[[x]],2, sum), "/")
                            res[is.na(res)] <- 0 
                            return(res)})      
    names(QS) <- stnms

    fl    <- fleets[[flnm]]
    
    sts   <- catchNames(fl)
    
    # The effort is restricted only by the stocks in 'stocks.restr'    
    # Remove the NA-s if any
    # if(any(fleets.ctrl[[flnm]][['stocks.restr']])){
    #   cat(paste("warning: there is at least one NA in  fleets.ctrl[['",flnm,"']][['stocks.restr']], and it has been removed, only the values different to NA will be maintained.\n", sep=""))
    #   
    #   fleets.ctrl[[flnm]][['stocks.restr']] <- fleets.ctrl[[flnm]][['stocks.restr']][!is.na(fleets.ctrl[[flnm]][['stocks.restr']])]
    # }
    # If the restrictors are missing => all the stocks restrict.
    if(is.null(fleets.ctrl[[flnm]][['stocks.restr']]) |  length(fleets.ctrl[[flnm]][['stocks.restr']]) == 0) {
      fleets.ctrl[[flnm]][['stocks.restr']] <- catchNames(fleets[[flnm]])
    }  
      
    sts <- fleets.ctrl[[flnm]][['stocks.restr']]
    
    
    mtnms <- names(fl@metiers)
    
    # flinfo: matrix with information on which metier catch which stock.
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')


    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    vc.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@vcost[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    effs <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))
    Cr.f <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))

    q.m <- alpha.m <- beta.m  <- ret.m <- wd.m <- wl.m <- pr.m <- vector('list', length(sts))
    names(q.m) <- names(alpha.m) <- names(beta.m) <- names(ret.m) <- names(wl.m) <- names(wd.m) <- names(pr.m) <-sts

    for(st in sts){     # q.m, alpha.m.... by metier but stock specific

        # identify the first metier that catch stock st
        mtst <- flinfo[[st]][2]
            
        age.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[1]]
        age.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[1]]
        age.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[1]]

        unit.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[3]]
        unit.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[3]]
        unit.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[3]]

        q.m[[st]]     <- array(0, dim = c(length(mtnms), length(age.q), length(unit.q),it),     dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1:it))
        alpha.m[[st]] <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), it), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1:it))
        beta.m[[st]]  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        ret.m[[st]]   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        wl.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        wd.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        pr.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        

        for(mt in mtnms){

            if(!(st %in% names(fl@metiers[[mt]]@catches))) next
                    
            q.m[[st]][mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE] 
            alpha.m[[st]][mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE] 
            beta.m[[st]][mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE] 
            ret.m[[st]][mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, drop = TRUE] 
            wl.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, drop = TRUE]
            wd.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, drop = TRUE]
            pr.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@price[,yr,,ss, drop = TRUE]
        }    
        
        Cr.f[st,] <- TAC[st,]*QS[[st]][flnm,]
        
        if (length(fleets.ctrl[[flnm]]$LandObl)==1){
          LO<-fleets.ctrl[[flnm]]$LandObl
        }else{
          LO<-fleets.ctrl[[flnm]]$LandObl[yr]
        }
        
        if(LO){
          if(fleets.ctrl[[flnm]]$LandObl_yearTransfer[yr-1] == TRUE){ # If landing Obligation = TRUE discount possible quota transfer from previous year.
            Cr.f[st,] <- Cr.f[st,] - fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[st,yr-1,]
            Cr.f[st,] <- ifelse(Cr.f[st,]<0, 0, Cr.f[st,])  # if lower than 0 , set it to 0.
          }
        }
        
    }
    
    
  #  browser()
            
    ## Update the effort-share using the defined model
    effortShare.fun <- fleets.ctrl[[flnm]][['effshare.model']]
    efs.m <- eval(call(effortShare.fun, Cr = Cr.f,  N = N, B = B, q.m = q.m, rho = rho, efs.m = efs.m, alpha.m, 
                       beta.m = beta.m, ret.m = ret.m, wl.m = wl.m, wd.m = wd.m, pr.m = pr.m, vc.m = vc.m,
                       season = ss, year = yr, fleet = fl, fleet.ctrl = fleets.ctrl[[flnm]], restriction = restriction))
    
    cat('Effort share: ', efs.m, ', sum:', apply(efs.m,2,sum), '\n')
    # Update the fleets object with the new effort share
    for(mt in names(fl@metiers))  fl@metiers[[mt]]@effshare[,yr,,ss] <-  efs.m[mt,] 
    
    for(st in sts){
        
        effort.fun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'effort', sep = '.')
        for(i in 1:it){
          
           if(!is.null(dim(rho))) rhoi <- rho[st,i,drop=F]
           else rhoi <- rho[st]
           
           Nst  <- array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
            
           effs[st, i] <-  eval(call(effort.fun, Cr = Cr.f[st,i],  N = Nst[,,i,drop=F], q.m = q.m[[st]][,,,i,drop=F], rho = rhoi,
                                efs.m = efs.m[,i,drop=F], alpha.m = alpha.m[[st]][,,,i,drop=F], beta.m = beta.m[[st]][,,,i,drop=F],
                                ret.m = ret.m[[st]][,,,i,drop=F], wl.m = wl.m[[st]][,,,i,drop=F], wd.m = wd.m[[st]][,,,i,drop=F],
                                restriction = restriction))
        }
    }
    
    if(LO == FALSE){
        # Choose the effort.
        if(length(fleets.ctrl[[flnm]]$effort.restr)==1){
          rule=fleets.ctrl[[flnm]]$effort.restr
        }else{
          rule=fleets.ctrl[[flnm]]$effort.restr[yr]
        }
        eff <- effRule.SMFB(effs = effs, prev.eff = matrix(fl@effort[,yr-1,,ss,drop=T],1,it), rule = rule)
        # Capacity restrictions.  
        eff <- capacityRest.SMFB(eff, c(fl@capacity[,yr,,ss,drop=T]))                                   
#         fleets[[flnm]]@effort[,yr,,ss] <- eff 
        fl@effort[,yr,,ss] <- eff 
    }
    else{ # landObl == TRUE
      eff <- numeric(it)
      discount_yrtransfer <- matrix(0,nst,it, dimnames = list(sts,1:it))
      ret.m.new <- ret.m # retention may change derived from minimis exemption.
      min_ctrl <- rep(FALSE, length(sts))
      names(min_ctrl) <- sts
      
      # Identify the stocks that are unable to 'receive' any extra TAC from others due to overfishing.
      stks_OF <- overfishing(biols, fleets, advice.ctrl, yr) # matrix[nst,it]
      
      
      # Identify the minimum effort and compare with capactity, if > capacity => eff = capacity and the algorithm finish.
        for(i in 1:it){
            Emin <- min(effs[,i])
            if(Emin > c(fl@capacity[,yr,,ss,,i,drop=T])){ 
              fl@effort[,yr,,ss,,i] <- fl@capacity[,yr,,ss,,i,drop=T] 
                next
            }
            else{ # Minimis, Quota transfer btw years and QuotaSwap.
              
              minimis <- fleets.ctrl[[flnm]]$LandObl_minimis # logical(ny)
              yrtrans <- fleets.ctrl[[flnm]]$LandObl_yearTransfer # logical(ny)
              
              Ni         <- lapply(N, function(x) array(x[,,,,,i], dim= c(dim(x)[c(1,3)],1)))
              q.m.i      <- lapply(q.m, function(x) x[,,,i,drop=F])
              alpha.m.i  <- lapply(alpha.m, function(x) x[,,,i,drop=F])
              beta.m.i   <- lapply(beta.m, function(x) x[,,,i,drop=F])
              wl.m.i     <- lapply(wl.m, function(x) x[,,,i,drop=F])
              wd.m.i     <- lapply(wd.m, function(x) x[,,,i,drop=F])
              ret.m.i    <- lapply(ret.m, function(x) x[,,,i,drop=F])
              K <- c(fl@capacity[,yr,,ss,,i,drop=T])
              
              if(!is.null(dim(rho))) rhoi <- rho[st,i]
              else rhoi <- rho[st]
              
              names(Ni) <- names(N)
              names(q.m.i) <- names(q.m.i) <- names(q.m.i) <- names(q.m.i) <- names(q.m.i) <- names(q.m.i) <- names(q.m)
              
             Cr.f_min_qt <- Cr.f
             eff_min_qt <- effs[, i]
              # Minimis and Quota transfer.
              if(minimis[yr] == TRUE | yrtrans[yr] == TRUE){
                              
                eff_min_qt <- numeric(length(Ni))
                names(eff_min_qt) <- stnms
                
                Cr.f_min_qt <- Cr.f
              
                for(st in sts){
         #         browser()
                  effort.fun <- ifelse(dim(Ni[[st]])[1] == 1, 'CobbDouglasBio.effort', 'CobbDouglasAge.effort')
                  # To calculate the final quota, the year transfer % needs to be applied to the original quota before
                  # discounting the quota used the pevious year and then discount this quota.
                  min_p <- fleets.ctrl[[flnm]]$LandObl_minimis_p[st,yr] # matrix(st,ny)
                  yrt_p <- fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[st,yr] # matrix(st,ny)
                  
                  
                  Cr.f_min_qt[st,i] <- (Cr.f[st,i] + fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[st,yr-1,i])*(1+min_p+yrt_p) - # The quota restriction is enhanced in the proportion allowed by minimis and year transfer.
                                        fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[st,yr-1,i]
                  
                  eff_min_qt[st] <-  eval(call(effort.fun, Cr = Cr.f_min_qt[st,i],  N = Ni[[st]], q.m = q.m.i[[st]], rho =rhoi,
                                       efs.m = efs.m[,i,drop=F], alpha.m = alpha.m.i[[st]], beta.m = beta.m.i[[st]],
                                        ret.m = ret.m.i[[st]], wl.m = wl.m.i[[st]], wd.m = wd.m.i[[st]],
                                        restriction = restriction)) # the restriction in landing obligation should be catch
                }
              }
              E1 <- min(eff_min_qt) # The effort resulting from minimis and year quota transfer examptions.
                                      # We will use this effort later to divide the extra catch, in discards (from minimis), year quota transfer 
                                      # to discount in the following year and quota swap (in this order)
              
              
              # Quota Swap
              if(!is.null(dim(rho))) rhoi <- rho[,i]
              else rhoi <- rho
              
              fcube_lo <- QuotaSwap(stknms = sts, E1, Cr.f = Cr.f[,i], Cr.f_exemp = Cr.f_min_qt[,i], N = Ni, B = B[,i,drop=F], efs.m = efs.m[,i,drop=F], q.m = q.m.i, alpha.m = alpha.m.i, beta.m = beta.m.i, 
                                      wl.m = wl.m.i, wd.m = wd.m.i, ret.m = ret.m.i, K = K, rho = rhoi, flnm = flnm, fleets.ctrl = fleets.ctrl, stks_OF = stks_OF[,i],approach = 'fcube')
      
              eff[i] <- fcube_lo$E
              fl@effort[,yr,,ss,,i] <- fcube_lo$E
                  cat('Effort after Landing Obligation Exemptions: ',fcube_lo$E, '\n')
              
              # Divide the extra catch, in discards (from minimis, only those derived from MLS), year quota transfer 
              # to discount in the following year and quota swap (in this order)
              # discount_yrtransfer must be discounted from the quota next year.
     
              catch_Elo <- fcube_lo$catch
              diff      <- catch_Elo[sts]/Cr.f[sts,i] #[nst]
              diff <- ifelse(Cr.f[sts,i]  == 0 & catch_Elo[sts] == 0, 0, diff)
              discount_yrtransfer[,i] <- ifelse(diff < 1 + fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr], 0, 
                                        ifelse((diff - fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr] - 1) < fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[,yr], 
                                               (diff - fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr] - 1),
                                                fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[,yr]))*Cr.f[,i]
              
              # update ret.m to account for the discards due to minimise exemption.
              for(st in sts){
             # if(flnm == 'MON_OT' & yr == 41)
              #  browser()
                # if discards due to size are higher than discards allowed by minimise, ret.m.i is not changed,
                # otherwise it is increased so that the total discards equal to min_p*Cr.f  
                
                Cr.f[st,i] <- ifelse(Cr.f[st,i] == 0, 1e-6, Cr.f[st,i])
                min_p <- fleets.ctrl[[flnm]]$LandObl_minimis_p[st,yr] # matrix(st,ny)
                yrt_p <- fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[st,yr] # matrix(st,ny)
                Ca <- fcube_lo$Ca[[st]] # catch at age in weight
                Da <- fcube_lo$Da[[st]]
                Ds <- sum(Da)                
                ret.m.new[[st]][,,,i] <- ret.m[[st]][,,,i] - ifelse(Ds/Cr.f[st,i] > min_p, 0, min_p- Ds/Cr.f[st,i])
                min_ctrl[st] <- ifelse(Ds/Cr.f[st,i]  > min_p, FALSE, TRUE)
              }
              
            }
        }
    
      # Update the retention curve according to minimis.
      if(any(min_ctrl)){
          sts_min <- names(which(min_ctrl))
     #     browser()
          for(mt in names(fl@metiers)){
           if(any(sts_min %in% catchNames(fl@metiers[[mt]]))){
             for(st in sts_min[which(sts_min %in% catchNames(fl@metiers[[mt]]))]){
        
  
               fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,] <- ret.m.new[[st]][mt,,,]
               fl@metiers[[mt]]@catches[[st]]@discards.sel[,yr,] <- 1-ret.m.new[[st]][mt,,,]
     
             }
           }    
          }
     
          fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[,yr,] <- discount_yrtransfer
      }
    }

          
    #    save(advice,alpha.m,B,beta.m,Cr.f,rho,eff,effs,efs.m,fleets.ctrl, 
    #         q.m,QS,QS.ss,TAC,TAC.yr, file = paste(flnm, file = '.RData', sep = ""))
   
   # Update the quota share of this step and the next one if the 
   # quota share does not coincide with the actual catch. (update next one only if s < ns).
   for(st in sts){

        yr.share       <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [it]
        ss.share       <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, it))# [it,ns]
        quota.share.OR <- matrix(t(yr.share*ss.share), ns, it)
        # The catch.
        catchFun <- fleets.ctrl[[flnm]][[st]][['catch.model']]
       Nst  <-  array(N[[st]][drop=T],dim = dim(N[[st]])[c(1,3,6)])
        catchD <- eval(call(catchFun, N = Nst,  E = eff, efs.m = efs.m, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], wd.m = wd.m[[st]], wl.m = wl.m[[st]], ret.m = ret.m[[st]]))
        itD <- ifelse(is.null(dim(catchD)), 1, length(dim(catchD)))
        catch <- apply(catchD, itD, sum)  # sum catch along all dimensions except iterations.
            
        quota.share    <- updateQS.SMFB(QS = quota.share.OR, TAC = TAC.yr[st,], catch = catch, season = ss)        # [ns,it]
                              
        fleets.ctrl$seasonal.share[[st]][flnm,yr,,] <- t(t(quota.share)/apply(quota.share, 2,sum)) #[ns,it], doble 't' to perform correctly de division between matrix and vector.
         
   }
  
  fleets[[flnm]] <- fl
    
    return(list(fleets = fleets, fleets.ctrl = fleets.ctrl))
}



#-------------------------------------------------
## GRAVITY MODEL TO UPDATE THE EFFORT SHARE
#-------------------------------------------------

gravity.flbeia <- function(Cr,  N, B, q.m, rho, efs.m, alpha.m, beta.m,
                    ret.m, wl.m, wd.m, pr.m, vc.m, season, year, fleet, fleet.ctrl, restriction = restriction,...){ 
  
  N0 <- lapply(names(N), function(x) array(N[[x]], dim = dim(N[[x]])[c(1,3,6)]))
  names(N0) <- names(N)
  
  if(fleet.ctrl$gravity.model == 'revenue'){  
    V.m  <- Reduce('+', lapply(names(q.m), function(x) 
                              apply(q.m[[x]]*(sweep(wl.m[[x]], 2:4, N0[[x]], "*")^beta.m[[x]])*ret.m[[x]]*pr.m[[x]],c(1,4),sum)))
    TotV <- apply(V.m,2,sum)
    res <- sweep(V.m, 2, TotV, "/")
  }else{
    if(fleet.ctrl$gravity.model == 'profit'){
        V.m  <- Reduce('+', lapply(names(q.m), function(x) 
            apply(q.m[[x]]*(sweep(wl.m[[x]], 2:4, N0[[x]], "*")^beta.m[[x]])*ret.m[[x]]*pr.m[[x]],c(1,4),sum)))
        TotV <- apply(V.m - vc.m,2,sum) 
        
        res <- sweep(V.m, 2, TotV, "/")
      
    }
    else stop('gravity.model argument must be equal to "profit" or "revenue')
  }
  
  trad <- ifelse(is.null(fleet.ctrl$gravity.tradition), 0, fleet.ctrl$gravity.tradition)
  
  res <- efs.m*trad+ res*(1-trad)
  

  
  return(res)
}



#-------------------------------------------------
## mlogit MODEL TO UPDATE THE EFFORT SHARE
#-------------------------------------------------
mlogit.flbeia <- function(Cr, N, B, q.m, rho, efs.m, alpha.m, 
                          beta.m, ret.m, wl.m, wd.m, pr.m, vc.m,
                          season, year, fleet, fleet.ctrl, restriction,...){
  
  
  ## step 1 
  predict.df <- make_RUM_predict_df(model = fleet.ctrl[['mlogit.model']], fleet = fleet, s = season)
  
  res <- efs.m
  res[] <- NA
  
  for(i in 1:dim(N[[1]])[6]){
  ## step 2 
    
    Ni         <- lapply(N, function(x) x[,,,,,i, drop=F])
    q.m.i      <- lapply(q.m, function(x) x[,,,i,drop=F])
    alpha.m.i  <- lapply(alpha.m, function(x) x[,,,i,drop=F])
    beta.m.i   <- lapply(beta.m, function(x) x[,,,i,drop=F])
    wl.m.i     <- lapply(wl.m, function(x) x[,,,i,drop=F])
    wd.m.i     <- lapply(wd.m, function(x) x[,,,i,drop=F])
    ret.m.i    <- lapply(ret.m, function(x) x[,,,i,drop=F])
    pr.m.i     <- lapply(pr.m, function(x) x[,,,i,drop=F])
    
    updated.df <- update_RUM_params(model = fleet.ctrl[['mlogit.model']], predict.df = predict.df, 
                                  fleet = fleet, covars = covars, season = season, year = year,
                                  N = Ni, q.m = q.m.i, wl.m = wl.m.i, beta.m = beta.m.i, ret.m = ret.m.i, pr.m = pr.m.i) 
    ## step 3 
    res[,i] <- predict_RUM(model = fleet.ctrl[['mlogit.model']], updated.df = updated.df)
  }
  

  return(res)
}


# ** make_RUM_predict_df **:  this makes the correctly formated dataframe over which to predict 
#     the effort shares. It requires the mlogit model, fleet object and season as input.
make_RUM_predict_df <- function(model = NULL, fleet = NULL, season) {
  
  ## Pass mlogit model object
  ## Pass fleet object
  
  mod.coefs <- names(coef(model)) ## Model coefficients
  
  ## 1. season - note, just return the season for which we're predicting
  seas <- if(any(grepl("season", mod.coefs))) { season } else { NA }
  
  ## 2. catch or catch rates
  C <- if(any(sapply(catchNames(fleet), grepl, mod.coefs))) {
    
    ## Return the catchnames that are in the coefficients
    catchNames(fleet)[unlist(sapply(catchNames(fleet), function(n) { any(grepl(n, mod.coefs))}))]
    
  } else { NA }
  
  ## 3. vcost 
  v    <- if(any(grepl("vcost", mod.coefs))) { -1 } else { NA }
  
  ## 4. effshare 
  e    <- if(any(grepl("effshare", mod.coefs))) { -1 } else { NA }
  
  ## Construct the dataframe
  predict.df <- expand.grid(metier = fleet@metiers@names, 
                            choice = "yes", 
                            season = as.numeric(seas), 
                            vcost = v, 
                            effshare = e,
                            stringsAsFactors = FALSE)
  ## Remove any columns with NAs, indicating variable not used
  predict.df <- predict.df[,which(sapply(predict.df, function(x) all(!is.na(x))))]
  
  ## Combine with the catch rate columns
  if(!all(is.na(C))) {	
    C.df <- as.data.frame(matrix(-1, ncol = length(C), nrow = nrow(predict.df)))
    colnames(C.df) <- C
    
    predict.df <- cbind(predict.df, C.df)
  }
  
  predict.df$index <- seq_len(nrow(predict.df)) 
  ## Use mFormula to define model form
  LD.predict <- mlogit:::mlogit.data(predict.df, choice = "choice", shape = "long",
                            alt.var = "metier", chid.var = "index")
  
  return(LD.predict)
}

# ** update_RUM_params **: For this I have tried to keep the inputs the same as for the gravity model. 
#                 Here, we update the data in the predict_df (from 1) with the values to predict over.
update_RUM_params <- function(model = NULL, predict.df, fleet, covars, season, year,
                              N, q.m, wl.m, beta.m, ret.m, pr.m) {
  
  ## Update the values in the predict.df
  
  ## 2. catch / catch rates - on same scale.
  ## Note, these should be updated based on the biomass increases, so we do a
  ## similar calculation as for the gravity model
  ## Here have to be careful as not all metiers may catch all stocks...
  
  if(any(sapply(catchNames(fleet), grepl, names(coef(model))))) {
    
    N0 <- lapply(names(N), function(x) array(N[[x]], dim = dim(N[[x]])[c(1,3,6)]))
    names(N0) <- names(N)
    
    ## This should be the catch rate per stock per metier ??
    CR.m   <- lapply(names(q.m), function(x) 
      cbind(stock = x,
            as.data.frame(
              apply(q.m[[x]]*(sweep(wl.m[[x]], 2:4, N0[[x]], "*")^beta.m[[x]])*ret.m[[x]]*pr.m[[x]],c(1,4),sum)
            )
      )
    )
    
    CR <- do.call(rbind, CR.m)
    
    for(st in unique(CR$stock)) {
      predict.df[,st] <- CR[CR$stock == st,2] 
    }
    predict.df[is.na(predict.df)] <- 0
    
  }
  
  # 3. vcost
  if("vcost" %in% colnames(predict.df)) {
    v <- do.call(rbind, lapply(fl@metiers, function(x) cbind(metier = x@name,as.data.frame(x@vcost[,year,,season]))))
    predict.df$vcost <- v$data
  }
  
  # 4. effort share - past effort share, y-1
  if("effshare" %in% colnames(predict.df)) {
    e <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@effshare[,year-1,,season]))))
    predict.df$effshare <- e$data
  }
  
  return(predict.df)
  
}


# ** predict_RUM ** : this function does the predictions and returns the effort shares.
predict_RUM <- function(model, updated.df) {
  
  ## Extract the model matrix and parameter coefficients
  mod.mat <- model.matrix(model$formula, data = updated.df)
  beta <- as.matrix(coef(model))
  
  ## Check the model matrix and coefficients are ordered correctly
  if(any(!colnames(mod.mat) == rownames(beta))) {
    stop("Model matrix and coefficients are not the same")
  }
  
  ## linear predictor long
  eta_long <- mod.mat %*% beta
  
  ## linear predictor wide
  eta_wide <- matrix(eta_long, ncol = length(unique(updated.df$metier)), byrow = TRUE)
  names(eta_wide) <- updated.df$metier 
  
  ## convert to a probability
  p_hat <- exp(eta_wide) / rowSums(exp(eta_wide))
  colnames(p_hat) <- updated.df$metier 
  p_hat <- as.data.frame(t(p_hat))
  
#  cat('Effort Share mlogit: ', p_hat[,1], '\n')
  
  return(p_hat[,1])
  
}



#-------------------------------------------------
## MARKOV MODEL TO UPDATE THE EFFORT SHARE
#-------------------------------------------------

Markov.flbeia <- function(Cr, N, B, q.m, rho, efs.m, alpha.m, 
                          beta.m, ret.m, wl.m, wd.m, pr.m, vc.m,
                          season, year, fleet, fleet.ctrl, restriction,...){
  
  
  ## step 1 
  predict.df <- make_Markov_predict_df(model = fleet.ctrl[['Markov.model']], fleet = fleet, s = season)
  
  res <- efs.m
  res[] <- NA
  
  for(i in 1:dim(N[[1]])[6]){
    ## step 2 
    
    Ni         <- lapply(N, function(x) x[,,,,,i, drop=F])
    q.m.i      <- lapply(q.m, function(x) x[,,,i,drop=F])
    alpha.m.i  <- lapply(alpha.m, function(x) x[,,,i,drop=F])
    beta.m.i   <- lapply(beta.m, function(x) x[,,,i,drop=F])
    wl.m.i     <- lapply(wl.m, function(x) x[,,,i,drop=F])
    wd.m.i     <- lapply(wd.m, function(x) x[,,,i,drop=F])
    ret.m.i    <- lapply(ret.m, function(x) x[,,,i,drop=F])
    pr.m.i     <- lapply(pr.m, function(x) x[,,,i,drop=F])
    
    updated.df <- update_Markov_params(model = fleet.ctrl[['Markov.model']], predict.df = predict.df, 
                                    fleet = fleet, covars = covars, season = season, year = year,
                                    N = Ni, q.m = q.m.i, wl.m = wl.m.i, beta.m = beta.m.i, ret.m = ret.m.i, pr.m = pr.m.i) 
    ## step 3 
  #  browser()
    res[,i] <- predict_Markov(model = fleet.ctrl[['Markov.model']], updated.df = updated.df, fleet = fleet, season = season, year = year)
  }
  
  
  return(res)
}


make_Markov_predict_df <- function(model = NULL, fleet = NULL, season) {
  
  ## Pass multinom model object
  ## Pass fleet object
  
  mod.coefs <- model$coefnames ## Model coefficients
  
  ## 1. season - note, just return the season for which we're predicting
  seas <- if(any(grepl("season", mod.coefs))) { season } else { NA }
  
  ## 2. catch or catch rates
  C <- if(any(sapply(catchNames(fleet), grepl, mod.coefs))) {
    
    ## Return the catchnames that are in the coefficients
    catchNames(fleet)[unlist(sapply(catchNames(fleet), function(n) { any(grepl(n, mod.coefs))}))]
    
  } else { NA }
  
  ## 3. vcost 
  v    <- if(any(grepl("vcost", mod.coefs))) { -1 } else { NA }
  
  ## 4. effshare 
  e    <- if(any(grepl("effshare", mod.coefs))) { -1 } else { NA }
  
  ## Construct the dataframe
  ## Note, we need the state from which vessels are coming
  predict.df <- expand.grid(state.tminus1 = fleet@metiers@names,
                            season = as.numeric(seas), 
                            vcost = v, 
                            effshare = e,
                            stringsAsFactors = FALSE)
  ## Remove any columns with NAs, indicating variable not used
  predict.df <- predict.df[,which(sapply(predict.df, function(x) all(!is.na(x))))]
  
  ## Correct attributes for prediction data
  if(!is.na(seas)) { 
    if(attr(model$terms, "dataClasses")[["season"]] == "factor") {
      predict.df$season <- as.factor(predict.df$season)
    }
  }
  
  ## Combine with the catch rate columns
  if(!all(is.na(C))) {	
    C.df <- as.data.frame(matrix(-1, ncol = length(C), nrow = nrow(predict.df)))
    colnames(C.df) <- C
    
    predict.df <- cbind(predict.df, C.df)
  }
  
  return(predict.df)
}



update_Markov_params <- function(model = NULL, predict.df, fleet, covars, season, year,
                                 N, q.m, wl.m, beta.m, ret.m, pr.m) {
  
  ## Update the values in the predict.df
  
  ## 2. catch / catch rates - on same scale.
  ## Note, these should be updated based on the biomass increases, so we do a
  ## similar calculation as for the gravity model
  ## Here have to be careful as not all metiers may catch all stocks...
  
  if(any(sapply(catchNames(fleet), grepl, model$coefnames))) {
    
    N0 <- lapply(names(N), function(x) array(N[[x]], dim = dim(N[[x]])[c(1,3,6)]))
    names(N0) <- names(N)
    
    ## This should be the catch rate per stock per metier ??
    CR.m   <- lapply(names(q.m), function(x) 
      cbind(stock = x,
            as.data.frame(
              apply(q.m[[x]]*(sweep(wl.m[[x]], 2:4, N0[[x]], "*")^beta.m[[x]])*ret.m[[x]]*pr.m[[x]],c(1,4),sum)
            )
      )
    )
    
    CR <- do.call(rbind, CR.m)
    
    for(st in unique(CR$stock)) {
      predict.df[,st] <- CR[CR$stock == st,2]  ## This will repeat, to ensure we get for each metier combinations
    }
    predict.df[is.na(predict.df)] <- 0
    
  }
  
  # 3. vcost
  if("vcost" %in% colnames(predict.df)) {
    v <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@vcost[,year,,season]))))
    predict.df$vcost <- v$data
  }
  
  # 4. effort share - past effort share, y-1
  if("effshare" %in% colnames(predict.df)) {
    e <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@effshare[,year-1,,season]))))
    predict.df$effshare <- e$data
  }
  
  return(predict.df)
  
}



predict_Markov <- function(model, updated.df, fleet, season, year) {
  
  # Transition probs
# browser()
  p_hat <- cbind(updated.df[c("state.tminus1")], nnet:::predict.multinom(model, updated.df, type = "probs"))
  p_hat_mat <- as.matrix(p_hat[,2:ncol(p_hat)])
  
  # past effort
  
  # New year
  if(season == 1) {
    last.season <- dims(fleet)[["season"]]
    cur.eff <- as.matrix(sapply(fleet@metiers, function(x) x@effshare[,year-1, , last.season]))
  }
  
  # Same year
  if(season > 1) {
    cur.eff <- as.matrix(sapply(fleet@metiers, function(x) x@effshare[, year, , season-1]))
  }
  
  new.share <- apply(p_hat_mat, 2, function(x) x %*% cur.eff)
  
  if(round(sum(new.share),6) != 1) {stop("Error - effort share does not sum to 1")}
  
  return(new.share)
  
}



