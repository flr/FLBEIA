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
  
    # 'year' dimension.
    # Dimnsions and fl
    fl    <- fleets[[flnm]]
    
    # The effort is restricted only by the stocks in 'stocks.restr'    
    # If the restrictors are missing => all the stocks restrict.
    #-----------------------------------------------------------------------------------
    if(is.null(fleets.ctrl[[flnm]][['stocks.restr']]) |  length(fleets.ctrl[[flnm]][['stocks.restr']]) == 0) {
      fleets.ctrl[[flnm]][['stocks.restr']] <- catchNames(fleets[[flnm]])
    }  
    sts   <- intersect(fleets.ctrl[[flnm]][['stocks.restr']], catchNames(fl))
   
    stnms <- names(biols)
    mtnms <- names(fl@metiers)
    nmt   <- length(mtnms)
    nst   <- length(biols)
    ns    <- dim(biols[[1]]@n)[4] 
    dimnms <- dimnames(biols[[1]]@n) 
    nit <- dim(biols[[1]]@n)[6]
    
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
    
    
    # Advice season for each stock
    adv.ss <- setNames( rep(NA,nst), stnms)
    for (st in stnms) adv.ss[st] <- ifelse(is.null(advice.ctrl[[st]][["adv.season"]]), ns, advice.ctrl[[st]][["adv.season"]]) # [nst]
    
    
    # Transform the FLR objects into list of arrays in order to be able to work with non-FLR
    list2env(FLObjs2S3_fleetSTD(biols = biols, fleets = fleets, advice = advice, covars = covars, 
                                biols.ctrl = biols.ctrl, fleets.ctrl = fleets.ctrl, BDs=BDs, 
                                flnm = flnm, yr = yr, ss = ss, iters = 1:nit, adv.ss), environment())
    
    ## Update the effort-share using the defined model
    effortShare.fun <- fleets.ctrl[[flnm]][['effshare.model']]
    efs.m <- eval(call(effortShare.fun, Cr = Cr.f,  N = N, B = B, q.m = q.m, rho = rho, efs.m = efs.m, alpha.m, 
                       beta.m = beta.m, ret.m = ret.m, wl.m = wl.m, wd.m = wd.m, pr.m = pr.m, vc.m = vc.m,
                       season = ss, year = yr, fleet = fl, fleet.ctrl = fleets.ctrl[[flnm]], restriction = restriction, covars=covars))
    
    cat('Effort share: ', efs.m, ', sum:', apply(efs.m,2,sum), '\n')
    # Update the fleets object with the new effort share
    for(mt in names(fl@metiers))  fl@metiers[[mt]]@effshare[,yr,,ss] <-  efs.m[mt,] 
    
           
     for(st in sts){     # q.m, alpha.m.... by metier but stock specific

        effort.fun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'effort', sep = '.')
        for(i in 1:nit){
          
           if(!is.null(dim(rho))) rhoi <- rho[,i,drop=F]
           else rhoi <- matrix(rho, length(stnms), 1, dimnames = list(stnms, 1))
 
           # Extract the i-th element from the lists. 
            Ni       <- lapply(setNames(sts, sts), function(x) array(N[[x]][,,i,drop=T], dim = c(dim(N[[x]])[c(1,2)],1)))
            q.mi     <- lapply(setNames(sts, sts),   function(x) q.m[[x]][,,,i,drop=F])
            beta.mi  <- lapply(setNames(sts, sts),   function(x) beta.m[[x]][,,,i,drop=F])
            alpha.mi <- lapply(setNames(sts, sts),   function(x) alpha.m[[x]][,,,i,drop=F])
            ret.mi   <- lapply(setNames(sts, sts),   function(x) ret.m[[x]][,,,i,drop=F])
            wl.mi    <- lapply(setNames(sts, sts),   function(x) wl.m[[x]][,,,i,drop=F])
            wd.mi    <- lapply(setNames(sts, sts),   function(x) wd.m[[x]][,,,i,drop=F])
        
            Nyri_1   <- lapply(setNames(sts, sts), function(x) array(Nyr_1[[x]][,,i,drop=T], dim = c(dim(Nyr_1[[x]])[c(1,2)],1)))
            Cyri_1   <- lapply(setNames(sts, sts), function(x) array(Cyr_1[[x]][,,i,drop=T], dim = c(dim(Cyr_1[[x]])[c(1,2)],1)))
            Cfyri_1  <- lapply(setNames(sts, sts), function(x) array(Cfyr_1[[x]][,,i,drop=T], dim = c(dim(Cfyr_1[[x]])[c(1,2)],1)))
            Myri_1   <- lapply(setNames(sts, sts), function(x) array(Myr_1[[x]][,,i,drop=T], dim = c(dim(Myr_1[[x]])[c(1,2)],1)))
            Mi       <- lapply(setNames(sts, sts), function(x) array(M[[x]][,,i,drop=T], dim = c(dim(M[[x]])[c(1,2)],1)))
    
            effs[st, i] <-  eval(call(effort.fun, Cr = Cr.f[,i, drop=F],  N = Ni, q.m = q.mi, rho = rhoi, efs.m = efs.m[,i,drop=F], 
                                alpha.m = alpha.mi, beta.m = beta.mi, ret.m = ret.mi, wl.m = wl.mi, wd.m = wd.mi,stknm=st,
                                restriction = restriction,  QS.groups = fleets.ctrl[[flnm]][['QS.groups']],
                                tac=TAC[,i,drop=F], Cyr_1 = Cyri_1, Nyr_1 = Nyri_1, Myr_1 = Myri_1,  M = Mi, Cfyr_1 = Cfyri_1))
        }
    }
    
    if(LO == FALSE){
        # Choose the effort.
        if(length(fleets.ctrl[[flnm]]$effort.restr)==1){
          rule=fleets.ctrl[[flnm]]$effort.restr
        }else{
          rule=fleets.ctrl[[flnm]]$effort.restr[yr]
        }
        eff <- effRule.SMFB(effs = effs, prev.eff = matrix(fl@effort[,yr-1,,ss,drop=T],1,nit), rule = rule)
        # Capacity restrictions.  
        eff <- capacityRest.SMFB(eff, c(fl@capacity[,yr,,ss,drop=T]))                                   
        fl@effort[,yr,,ss] <- eff 
    }
    else{ # landObl == TRUE
      eff <- numeric(nit)
      discount_yrtransfer <- matrix(0,length(sts),nit, dimnames = list(sts,1:nit))
      ret.m.new <- ret.m # retention may change derived from minimis exemption.
      min_ctrl <- rep(FALSE, length(sts))
      names(min_ctrl) <- sts
      
      # Identify the stocks that are unable to 'receive' any extra TAC from others due to overfishing.
      stks_OF <- overfishing(biols, fleets, advice.ctrl, yr) # matrix[nst,nit]
      
      
      # Identify the minimum effort and compare with capactity, if > capacity => eff = capacity and the algorithm finish.
        for(i in 1:nit){
            Emin <- min(effs[,i])
            if(Emin > c(fl@capacity[,yr,,ss,,i,drop=T])){ 
              fl@effort[,yr,,ss,,i] <- fl@capacity[,yr,,ss,,i,drop=T] 
                next
            }
            else{ # Minimis, Quota transfer btw years and QuotaSwap.
              
              minimis <- fleets.ctrl[[flnm]]$LandObl_minimis # logical(ny)
              yrtrans <- fleets.ctrl[[flnm]]$LandObl_yearTransfer # logical(ny)
              
              
              if(!is.null(dim(rho))) rhoi <- rho[,i,drop=F]
              else rhoi <- matrix(rho, length(stnms), 1, dimnames = list(stnms, 1))
              
              # Extract the i-th element form the lists. 
              Ni       <- lapply(setNames(sts, sts), function(x) array(N[[x]][,,i,drop=T], dim = c(dim(N[[x]])[c(1,3)],1)))
              q.mi     <- lapply(setNames(sts, sts),   function(x) q.m[[x]][,,,i,drop=F])
              beta.mi  <- lapply(setNames(sts, sts),   function(x) beta.m[[x]][,,,i,drop=F])
              alpha.mi <- lapply(setNames(sts, sts),   function(x) alpha.m[[x]][,,,i,drop=F])
              ret.mi   <- lapply(setNames(sts, sts),   function(x) ret.m[[x]][,,,i,drop=F])
              wl.mi    <- lapply(setNames(sts, sts),   function(x) wl.m[[x]][,,,i,drop=F])
              wd.mi    <- lapply(setNames(sts, sts),   function(x) wd.m[[x]][,,,i,drop=F])
              
              K <- c(fl@capacity[,yr,,ss,,i,drop=T])
      
              Cr.f_min_qt <- Cr.f
              eff_min_qt <- effs[, i]
              # Minimis and Quota transfer.
              if(minimis[yr] == TRUE | yrtrans[yr] == TRUE){
                              
                eff_min_qt <- numeric(length(Ni))
                names(eff_min_qt) <- sts
                
                Cr.f_min_qt <- Cr.f
              
                for(st in sts){

                  if(!is.null(dim(rho))) rhoi <- rho[,i,drop=F]
                  else rhoi <- matrix(rho, length(stnms), 1, dimnames = list(stnms, 1))
                  
                  # Extract the i-th element form the lists. 
                  Ni       <- lapply(setNames(sts, sts), function(x) array(N[[x]][,,i,drop=T], dim = c(dim(N[[x]])[c(1,3)],1)))
                  q.mi     <- lapply(setNames(sts, sts),   function(x) q.m[[x]][,,,i,drop=F])
                  beta.mi  <- lapply(setNames(sts, sts),   function(x) beta.m[[x]][,,,i,drop=F])
                  alpha.mi <- lapply(setNames(sts, sts),   function(x) alpha.m[[x]][,,,i,drop=F])
                  ret.mi   <- lapply(setNames(sts, sts),   function(x) ret.m[[x]][,,,i,drop=F])
                  wl.mi    <- lapply(setNames(sts, sts),   function(x) wl.m[[x]][,,,i,drop=F])
                  wd.mi    <- lapply(setNames(sts, sts),   function(x) wd.m[[x]][,,,i,drop=F])
                  
                  Nyri_1   <- lapply(setNames(sts, sts), function(x) array(Nyr_1[[x]][,,i,drop=T], dim = c(dim(Nyr_1[[x]])[c(1,2)],1)))
                  Cyri_1   <- lapply(setNames(sts, sts), function(x) array(Cyr_1[[x]][,,i,drop=T], dim = c(dim(Cyr_1[[x]])[c(1,2)],1)))
                  Cfyri_1  <- lapply(setNames(sts, sts), function(x) array(Cfyr_1[[x]][,,i,drop=T], dim = c(dim(Cfyr_1[[x]])[c(1,2)],1)))
                  Myri_1   <- lapply(setNames(sts, sts), function(x) array(Myr_1[[x]][,,i,drop=T], dim = c(dim(Myr_1[[x]])[c(1,2)],1)))
                  Mi       <- lapply(setNames(sts, sts), function(x) array(M[[x]][,,i,drop=T], dim = c(dim(M[[x]])[c(1,2)],1)))
                 
                  effort.fun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'effort', sep = '.')
                  # To calculate the final quota, the year transfer % needs to be applied to the original quota before
                  # discounting the quota used the pevious year and then discount this quota.
                  min_p <- fleets.ctrl[[flnm]]$LandObl_minimis_p[st,yr] # matrix(st,ny)
                  yrt_p <- fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[st,yr] # matrix(st,ny)
                  
                  
                  Cr.f_min_qt[st,i] <- (Cr.f[st,i] + fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[st,yr-1,i])*(1+min_p+yrt_p) - # The quota restriction is enhanced in the proportion allowed by minimis and year transfer.
                                        fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[st,yr-1,i]
                  
                  eff_min_qt[st] <-  eval(call(effort.fun, Cr = Cr.f[,i, drop=F],  N = Ni, q.m = q.mi, rho = rhoi, efs.m = efs.m[,i,drop=F], 
                                               alpha.m = alpha.mi, beta.m = beta.mi, ret.m = ret.mi, wl.m = wl.mi, wd.m = wd.mi,stknm=st,
                                               restriction = restriction,  QS.groups = fleets.ctrl[[flnm]][['QS.groups']],
                                               tac=TAC[,i,drop=F], Cyr_1 = Cyri_1, Nyr_1 = Nyri_1, Myr_1 = Myri_1,  M = Mi, Cfyr_1 = Cfyri_1))
                  
            
                }
              }
              E1 <- min(eff_min_qt) # The effort resulting from minimis and year quota transfer examptions.
                                      # We will use this effort later to divide the extra catch, in discards (from minimis), year quota transfer 
                                      # to discount in the following year and quota swap (in this order)
              
              
              # Quota Swap
              if(!is.null(dim(rho))) rhoi <- rho[,i,drop=F]
              else rhoi <- matrix(rho, length(stnms), 1, dimnames = list(stnms, 1))
              
              fcube_lo <- QuotaSwap(stknms = sts, E1, Cr.f = Cr.f[,i], Cr.f_exemp = Cr.f_min_qt[,i], N = Ni, B = B[,i,drop=F], efs.m = efs.m[,i,drop=F], q.m = q.mi, alpha.m = alpha.mi, beta.m = beta.mi, 
                                      wl.m = wl.mi, wd.m = wd.mi, ret.m = ret.mi, K = K, rho = rhoi, flnm = flnm, fleets.ctrl = fleets.ctrl, stks_OF = stks_OF[,i],approach = 'fcube')
      
              eff[i] <- fcube_lo$E
              fl@effort[,yr,,ss,,i] <- fcube_lo$E
                  cat('Effort after Landing Obligation Exemptions: ',fcube_lo$E, '\n')
              
              # Divide the extra catch, in discards (from minimis, only those derived from MLS), year quota transfer 
              # to discount in the following year and quota swap (in this order)
              # discount_yrtransfer must be discounted from the quota next year.
     
              catch_Elo <- fcube_lo$catch
              diff      <- catch_Elo[sts]/Cr.f[sts,i] #[nst]
              diff <- ifelse(Cr.f[sts,i]  == 0 & catch_Elo[sts] == 0, 0, diff)
              discount_yrtransfer[,i] <- ifelse(diff < 1 + fleets.ctrl[[flnm]]$LandObl_minimis_p[sts,yr], 0, 
                                        ifelse((diff - fleets.ctrl[[flnm]]$LandObl_minimis_p[sts,yr] - 1) < fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[sts,yr], 
                                               (diff - fleets.ctrl[[flnm]]$LandObl_minimis_p[sts,yr] - 1),
                                                fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[sts,yr]))*Cr.f[,i]
              
              # update ret.m to account for the discards due to minimise exemption.
              for(st in sts){
              # if discards due to size are higher than discards allowed by minimise, ret.m.i is not changed,
              # otherwise nit is increased so that the total discards equal to min_p*Cr.f  
                
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
     
          fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[sts,yr,] <- discount_yrtransfer
      }
    }

   
   # Update the quota share of this step and the next one if the 
   # quota share does not coincide with the actual catch. (update next one only if s < ns).
   for(st in sts){

        if (adv.ss[st] == ns) {
          yr.share <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [nit]
          ss.share <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, nit)) # [nit,ns]
        } else {
          ss1 <- (adv.ss[st]+1):ns
          ss2 <- 1:adv.ss[st]
          
          if (ss <= adv.ss[st]) {
            yr.share <- advice$quota.share[[st]][flnm,yr-1,, drop=T]      # [nit]
            ss.share <- cbind( t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,ss2, drop=T], length(ss2), nit)), 
                               t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr-1,,ss1, drop=T], length(ss1), nit))) # [nit,ns]
          } else {
            yr.share <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [nit]
            ss.share <- cbind( t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr+1,,ss2, drop=T], length(ss2), nit)), 
                               t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,ss1, drop=T], length(ss1), nit))) # [nit,ns]
          }
        }
     
        quota.share.OR <- matrix(t(yr.share*ss.share), ns, nit)
        
        # The catch.
        catchFun <- fleets.ctrl[[flnm]][[st]][['catch.model']]
        
        catchD <- array(NA, dim=dim(q.m[[st]]))
        
        for(i in 1:nit){
   
          if(is.null(dim(rho)))   rhoi <- rho
          if(length(dim(rho))==2) rho <- rho[st,i]
          Nyri_1   <- lapply(setNames(sts, sts), function(x) array(Nyr_1[[x]][,,i,drop=T], dim = c(dim(Nyr_1[[x]])[c(1,2)],1)))
          Cyri_1   <- lapply(setNames(sts, sts), function(x) array(Cyr_1[[x]][,,i,drop=T], dim = c(dim(Cyr_1[[x]])[c(1,2)],1)))
          Cfyri_1  <- lapply(setNames(sts, sts), function(x) array(Cfyr_1[[x]][,,i,drop=T], dim = c(dim(Cfyr_1[[x]])[c(1,2)],1)))
          Myri_1   <- lapply(setNames(sts, sts), function(x) array(Myr_1[[x]][,,i,drop=T], dim = c(dim(Myr_1[[x]])[c(1,2)],1)))
          Mi       <- lapply(setNames(sts, sts), function(x) array(M[[x]][,,i,drop=T], dim = c(dim(M[[x]])[c(1,2)],1)))
          
          #browser()
          
           catchD[,,,i] <- eval(call(catchFun, Cr=Cr.f[st,i],N = Ni[[st]],  E = eff[i], efs.m = efs.m[,i,drop=FALSE], q.m = q.m[[st]][,,,i,drop=FALSE], 
                            alpha.m = alpha.m[[st]][,,,i,drop=FALSE], beta.m = beta.m[[st]][,,,i,drop=FALSE], wd.m = wd.m[[st]][,,,i,drop=FALSE],
                            wl.m = wl.m[[st]][,,,i,drop=FALSE], ret.m = ret.m[[st]][,,,i,drop=FALSE], rho = rho,
                            tac=TAC[st,i], Cyr_1 = Cyri_1[[st]], Nyr_1 = Nyri_1[[st]], Myr_1 = Myri_1[[st]],  M = Mi[[st]], 
                            Cfyr_1 = Cfyri_1[[st]]))
         }
        itD <- ifelse(is.null(dim(catchD)), 1, length(dim(catchD)))
        catch <- apply(catchD, itD, sum)  # sum catch along all dimensions except iterations.
            
        quota.share     <- updateQS.SMFB(QS = quota.share.OR, TAC = TAC.yr[st,], catch = catch, season = ss, adv.season = adv.ss[st]) # [ns,nit]
        
        quota.share.NEW <- t(t(quota.share)/apply(quota.share, 2,sum)) #[ns,nit] double 't' to perform correctly the division between matrix and vector. 
        
        if (adv.ss[st] == ns) {
          fleets.ctrl$seasonal.share[[st]][flnm,yr,,] <- quota.share.NEW
        } else {
          if (ss <= adv.ss[st]) {
            fleets.ctrl$seasonal.share[[st]][flnm,yr-1,,ss1,] <- quota.share.NEW[ss1,]
            fleets.ctrl$seasonal.share[[st]][flnm,yr,,ss2,]   <- quota.share.NEW[ss2,]
          } else {
            fleets.ctrl$seasonal.share[[st]][flnm,yr,,ss1,]   <- quota.share.NEW[ss1,]
            fleets.ctrl$seasonal.share[[st]][flnm,yr+1,,ss2,] <- quota.share.NEW[ss2,]
          }
        }
         
   }
  
  fleets[[flnm]] <- fl
    
    return(list(fleets = fleets, fleets.ctrl = fleets.ctrl))
}

#-------------------------------------------------
## GRAVITY MODEL TO UPDATE THE EFFORT SHARE
#-------------------------------------------------

gravity.flbeia <- function(Cr,  N, B, q.m, rho, efs.m, alpha.m, beta.m,
                    ret.m, wl.m, wd.m, pr.m, vc.m, season, year, fleet, fleet.ctrl, restriction = restriction,...){ 
  
  N0 <- N
  
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
  predict.df <- make_RUM_predict_df(model = fleet.ctrl[['mlogit.model']], fleet = fleet, season = season)
  
  res <- efs.m
  res[] <- NA
  
  for(i in 1:dim(N[[1]])[3]){
  ## step 2 
    
    Ni         <- lapply(N, function(x) x[,,i, drop=F])
    q.m.i      <- lapply(q.m, function(x) x[,,,i,drop=F])
    alpha.m.i  <- lapply(alpha.m, function(x) x[,,,i,drop=F])
    beta.m.i   <- lapply(beta.m, function(x) x[,,,i,drop=F])
    wl.m.i     <- lapply(wl.m, function(x) x[,,,i,drop=F])
    wd.m.i     <- lapply(wd.m, function(x) x[,,,i,drop=F])
    ret.m.i    <- lapply(ret.m, function(x) x[,,,i,drop=F])
    pr.m.i     <- lapply(pr.m, function(x) x[,,,i,drop=F])
    
    updated.df <- update_RUM_params(model = fleet.ctrl[['mlogit.model']], predict.df = predict.df, 
                                  fleet = fleet, covars = covars, season = season, year = year,
                                  N = Ni, q.m = q.m.i, wl.m = wl.m.i, beta.m = beta.m.i, ret.m = ret.m.i, pr.m = pr.m.i,
				  iter = i) 
    ## step 3

	# If all of the catch.q for a given metier are zero, that metier is closed.
	# so to work out which metier are closed
	met.close <- apply(do.call(rbind, lapply(q.m.i, function(x) apply(x==0,1,all))),2,all)
	met.close <- ifelse(identical(names(which(met.close == TRUE)), character(0)), NA,
		     names(which(met.close == TRUE)))

     res[,i] <- predict_RUM(model = fleet.ctrl[['mlogit.model']], updated.df = updated.df, season, close = met.close)
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
  ## Determine if a factor or numeric
    if(!is.na(seas)) {
        
    if(any(class(model.frame(model)$season) == "numeric")) { 
        seas <- as.numeric(seas) } else { 
        seas <- as.factor(seas)
        }
        
    ## If season is a factor, we need to include the other seasons for contrast
    if(class(seas) == "factor") {
    seas <- as.factor(1:max(as.numeric(as.character(model.frame(model)$season)), na.rm = T))
    }
        }
 
  
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
                            choice = c(TRUE,FALSE), 
                            season = seas, 
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
  LD.predict <- mlogit::mlogit.data(predict.df, choice = "choice", shape = "long",
                            alt.var = "metier", chid.var = "index")
  
  return(LD.predict)
}

# ** update_RUM_params **: For this I have tried to keep the inputs the same as for the gravity model. 
#                 Here, we update the data in the predict_df (from 1) with the values to predict over.
update_RUM_params <- function(model = NULL, predict.df, fleet, covars, season, year,
                              N, q.m, wl.m, beta.m, ret.m, pr.m, iter) {
  
  ## Update the values in the predict.df
  
  ## 2. catch / catch rates - on same scale.
  ## Note, these should be updated based on the biomass increases, so we do a
  ## similar calculation as for the gravity model
  ## Here have to be careful as not all metiers may catch all stocks...
  
  if(any(sapply(catchNames(fleet), grepl, names(coef(model))))) {
    
    N0 <- N
    
    ## catch rate per stock per metier 
    CR.m   <- lapply(names(q.m), function(x) 
      cbind(stock = x,
            as.data.frame(
              apply(q.m[[x]]*(sweep(wl.m[[x]], 2:4, N0[[x]], "*")^beta.m[[x]])*ret.m[[x]],c(1,4),sum)
            )
      )
    )
    
    CR <- do.call(rbind, CR.m)
    
    for(st in unique(CR$stock)) {
      predict.df[,st] <- CR[CR$stock == st,2] 
    }
    predict.df[is.na(predict.df),] <- 0
    
  }
  
  # 3. vcost
  if("vcost" %in% colnames(predict.df)) {
    v <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@vcost[,year,,season,,iter]))))
    predict.df$vcost <- v$data
  }
  
  # 4. effort share - past effort share, y-1
  if("effshare" %in% colnames(predict.df)) {
    e <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@effshare[,year-1,,season,,iter]))))
    predict.df$effshare <- e$data
  }
  
  return(predict.df)
  
}


# ** predict_RUM ** : this function does the predictions and returns the effort shares.
predict_RUM <- function(model, updated.df, season, close) {
 

  ## Just the predictions we're interested in...
  updated.df <- updated.df[updated.df$choice == TRUE &
			  updated.df$season == season,]

  ## Extract the model matrix and parameter coefficients
  mod.mat <- model.matrix(mlogit::mFormula(model$formula), data = updated.df)
  beta <- as.matrix(coef(model))
  
  ## Check the model matrix and coefficients are ordered correctly
  if(any(!colnames(mod.mat) == rownames(beta))) {
    stop("Model matrix and coefficients are not the same")
  }
    
  ## If season is a factor, we want to exclude these options and just get the 
  ## predictions for the relevant season. Note if season is a numeric, the model 
  ## matrix already only includes the right season
    
  if(any(grepl("season", colnames(mod.mat)))) {
      
  if(any(class(model.frame(model)$season) == "factor")) { 
      seas <- 1:max(as.numeric(as.character(model.frame(model)$season)),na.rm=T)
      toRemove <- paste0("season", seas[!seas %in% season])
      
      # remove from mod.mat
      mod.mat <- mod.mat[,!colnames(mod.mat) %in% grep(paste(toRemove, collapse = "|"), colnames(mod.mat), value = T)]
      # remove from beta
      beta <- beta[!rownames(beta) %in% grep(paste(toRemove, collapse = "|"), rownames(beta), value = T),]
  }
      
  }
    
  ## linear predictor long
  eta_long <- mod.mat %*% beta
  
  ## linear predictor wide
  eta_wide <- matrix(eta_long, ncol = length(unique(updated.df$metier)), byrow = TRUE)
  names(eta_wide) <- updated.df$metier 
  
  ## Implement spatial closures
  eta_wide[names(eta_wide) %in% close] <- -Inf
  
  ## convert to a probability
  p_hat <- exp(eta_wide) / rowSums(exp(eta_wide))
  colnames(p_hat) <- unique(updated.df$metier)
  p_hat <- as.data.frame(t(p_hat))
  
  return(p_hat[,1])
  
}



#-------------------------------------------------
## MARKOV MODEL TO UPDATE THE EFFORT SHARE
#-------------------------------------------------

Markov.flbeia <- function(Cr, N, B, q.m, rho, efs.m, alpha.m, 
                          beta.m, ret.m, wl.m, wd.m, pr.m, vc.m,
                          season, year, fleet, fleet.ctrl, restriction,...){
  
  args   <- list(...)
  covars <- args$covars
  
  ## step 1 
  predict.df <- make_Markov_predict_df(model = fleet.ctrl[['Markov.model']], fleet = fleet, season = season)
  
  res <- efs.m
  res[] <- NA
  
  for(i in 1:dim(N[[1]])[3]){
    ## step 2 
    Ni         <- lapply(N, function(x) x[,,i, drop=F])
    q.m.i      <- lapply(q.m, function(x) x[,,,i,drop=F])
    alpha.m.i  <- lapply(alpha.m, function(x) x[,,,i,drop=F])
    beta.m.i   <- lapply(beta.m, function(x) x[,,,i,drop=F])
    wl.m.i     <- lapply(wl.m, function(x) x[,,,i,drop=F])
    wd.m.i     <- lapply(wd.m, function(x) x[,,,i,drop=F])
    ret.m.i    <- lapply(ret.m, function(x) x[,,,i,drop=F])
    pr.m.i     <- lapply(pr.m, function(x) x[,,,i,drop=F])
      
   
    updated.df <- update_Markov_params(model = fleet.ctrl[['Markov.model']], predict.df = predict.df, 
                                    fleet = fleet, covars = covars, season = season, year = year,
                                    N = Ni, q.m = q.m.i, wl.m = wl.m.i, beta.m = beta.m.i, ret.m = ret.m.i, pr.m = pr.m.i, iter = i) 
    ## step 3 
    
    # If all of the catch.q for a given metier are zero, that metier is closed.
    # so to work out which metier are closed
    met.close <- apply(do.call(rbind, lapply(q.m.i, function(x) apply(x==0,1,all))),2,all)
    met.close <- ifelse(identical(names(which(met.close == TRUE)), character(0)), NA,
		     names(which(met.close == TRUE)))

    res[,i] <- predict_Markov(model = fleet.ctrl[['Markov.model']], updated.df = updated.df, fleet = fleet, season = season, year = year, close = met.close, iter = i)
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
                            season = seas, 
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
                                 N, q.m, wl.m, beta.m, ret.m, pr.m, iter) {
  
  ## Update the values in the predict.df
  
  ## 2. catch / catch rates - on same scale.
  ## Note, these should NOT be updated based on the biomass increases, 
  ## we take these from the previous season (as it should of been fitted)
  
  if(any(sapply(catchNames(fleet), grepl, model$coefnames))) {
    
## old method	  
#    N0 <- N
    ## This should be the catch rate per stock per metier ??
#    CR.m   <- lapply(names(q.m), function(x) 
#      cbind(stock = x,
#            as.data.frame(
 #             apply(q.m[[x]]*(sweep(wl.m[[x]], 2:4, N0[[x]], "*")^beta.m[[x]])*ret.m[[x]],c(1,4),sum)
  #          )
  #    )
  #  )
    
 #   CR <- do.call(rbind, CR.m)

 ## New method with lagged data
# if first season, last season previous year
  year_lag <- ifelse(season == 1, year-1, year)   
  seas_lag <- ifelse(season == 1, dim(fleet@effort)[4], season-1)

## Get the landings in last season
land <- do.call(rbind, lapply(fleet@metiers, function(m) {
do.call(rbind, lapply(m@catches, function(x) cbind(metier = m@name, stock= x@name,as.data.frame(x@landings[,year_lag,,seas_lag, , iter]))))
}))

## Get the metier effort last season
 eff <- do.call(rbind, lapply(fleet@metiers, function(x) { 
 cbind(metier = x@name,as.data.frame(x@effshare[,year_lag,,seas_lag, , iter]))}))

eff$data <- as.data.frame(fleet@effort[,year_lag,,seas_lag,,iter])$data * eff$data

# combine the effort and landings and calculate the lpue
land$effort <- eff$data[match(land$metier, eff$metier)]
land$lpue   <- land$data / land$effort 


  for(st in colnames(predict.df)) {
	 
     predict.df[predict.df$state.tminus1 %in% land[land$stock == st, "metier"],st] <- land[land$stock == st, "lpue"]
     #predict.df[,st] <- land[land$stock == st, "lpue"]
     # CR[CR$stock == st,2]  ## This will repeat, to ensure we get for each metier combinations
    }
    predict.df[is.na(predict.df)] <- 0
  }
  
  # 3. vcost
  if("vcost" %in% colnames(predict.df)) {
    v <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@vcost[,year_lag,,seas_lag, , iter]))))
    predict.df$vcost <- v$data
  }
  
  # 4. effort share - past effort share, y-1
  if("effshare" %in% colnames(predict.df)) {
    e <- do.call(rbind, lapply(fleet@metiers, function(x) cbind(metier = x@name,as.data.frame(x@effshare[,year_lag,,seas_lag, , iter]))))
    predict.df$effshare <- e$data
  }
  
  return(predict.df)
  
}



predict_Markov <- function(model, updated.df, fleet, season, year, close, iter = i) {
  
  # Transition probs
  p_hat <- cbind(updated.df[c("state.tminus1")], nnet:::predict.multinom(model, updated.df, type = "probs"))
  p_hat_mat <- as.matrix(p_hat[,2:ncol(p_hat)])
  
  ## Implement spatial closures
  p_hat_mat[,colnames(p_hat_mat) %in% close] <- 0
  p_hat_mat <- p_hat_mat / rowSums(p_hat_mat, na.rm = TRUE)
  
  # past effort
  
  # New year
  if(season == 1) {
    last.season <- dims(fleet)[["season"]]
    cur.eff <- as.matrix(sapply(fleet@metiers, function(x) x@effshare[,year-1, , last.season,, iter]))
  }
  
  # Same year
  if(season > 1) {
    cur.eff <- as.matrix(sapply(fleet@metiers, function(x) x@effshare[, year, , season-1,, iter]))
  }
  
  new.share <- apply(p_hat_mat, 2, function(x) x %*% cur.eff)
  
  if(round(sum(new.share),6) != 1) {stop("Error - effort share does not sum to 1")}
  
  return(new.share)
  
}

                             
