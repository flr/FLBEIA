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
    
#if(flnm == 'GNS_FR') browser()
    
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
                                flnm = flnm, yr = yr, ss = ss, iters = 1:nit, adv.ss = adv.ss), environment())
    
#    if(flnm == 'GNS_FR') browser()
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
     
        # if advice previous to final season and no more year in the object --> do not update
        if (ss == ns & adv.ss[st] < ns & yr ==  dim(fleets.ctrl$seasonal.share[[st]])[2])
          next()
       
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
          Ni       <- lapply(setNames(sts, sts), function(x) array(N[[x]][,,i,drop=T], dim = c(dim(N[[x]])[c(1,2)],1)))
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

