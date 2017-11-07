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
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SMFB_LO(fleets, biols, covars, fleets.ctrl, year = 1, season = 1)
#-------------------------------------------------------------------------------
SMFB <- function(fleets, biols, BDs, covars, advice, biols.ctrl, fleets.ctrl, advice.ctrl, flnm, year = 1, season = 1,...){
    
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
     
#     if(flnm == "ANK_OT") browser()


    
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
    mtnms <- names(fl@metiers)
    
    # flinfo: matrix with information on which metier catch which stock.
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')


    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    effs <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))
    Cr.f <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))

    q.m <- alpha.m <- beta.m  <- ret.m <- wd.m <- wl.m <-vector('list', length(sts))
    names(q.m) <- names(alpha.m) <- names(beta.m) <- names(ret.m) <- names(wl.m) <- names(wd.m) <- sts

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


        for(mt in mtnms){

            if(!(st %in% names(fl@metiers[[mt]]@catches))) next
                    
            q.m[[st]][mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE] 
            alpha.m[[st]][mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE] 
            beta.m[[st]][mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE] 
            ret.m[[st]][mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, drop = TRUE] 
            wl.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, drop = TRUE]
            wd.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, drop = TRUE]
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
                # otherwise it is increases so that the total discards equal to min_p*Cr.f  
                
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

