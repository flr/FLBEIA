#-------------------------------------------------------------------------------
#        Functions to calculate and update  *effort*, *effshare*, *catch.q*..
#    i.e the fleets' parameters neccesary to calculate the landings and 
#    discards at age.
#
#    - 'fixedEffort' - (All Parameters Are Given) - so the function just returns the 
#        given FLFleetsExt object. It is used to maintain generallity.
#    - 'SMFB' - (Simple mixed fisheries behaviour). - Everything constant e
#        except effort that is updated based on landings or catch share. 
#       (multiple TACs so min, max Effort options are applied)
#    - 'SSFB' - (Simple sequential fisheries behaviour). - Everything constant 
#        except effort that is updated based on landings or catch share 
#       (condition: each metier of the fleet targets only one stock)
#
# Dorleta GarcYYYa
# Created: 28/10/2010 12:33:04
# Changed: 2011-02-28 16:17:37 (ssanchez)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# fixedEffort(fleets, biols, covars, fleets.ctrl, year = 1, season = 1)
#-------------------------------------------------------------------------------
fixedEffort <- function(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1){
    return(list(fleets = fleets, fleets.ctrl = fleets.ctrl))
}

#-------------------------------------------------------------------------------
# SMFB(fleets, biols, covars, fleets.ctrl, year = 1, season = 1)
#-------------------------------------------------------------------------------
SMFB <- function(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1){
    
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
    if(! all(sapply(names(fleets), function(x) fleets.ctrl[[x]]$restriction %in% c('catch', 'landings'))))
        stop("fleets.ctrl$restriction must be equal to 'catch' or 'landings'")
     
    # Dimensions.
    nst <- length(biols);          stnms <- names(biols)
    ns  <- dim(biols[[1]]@n)[4]
    it  <- dim(biols[[1]]@n)[6]
    flnms <- names(fleets)
    
    # Data
    B    <- matrix(t(sapply(stnms, function(x){   # biomass in the middle of the season  [nst,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss, drop=T])
                                else return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss, drop=T])})) , nst,it, dimnames = list(stnms, 1:it))

    Ba   <- lapply(stnms, function(x){   # biomass at age in the middle  of the season, list elements: [na,1,nu,1,1,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return((biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2))[,yr,,ss, drop = FALSE])
                                else return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss])})
    names(Ba) <- stnms
    

                      
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
    TAC.yr   <- matrix(advice$TAC[,yr,drop=T], nst,it, dimnames = list(stnms, 1:it))                       # [nst,it]
    CT       <- fleets.ctrl$catch.threshold[,yr,,ss, drop=T]  # [ns,it]
    QS.ss    <- matrix(t(sapply(stnms, function(x) apply(QS[[x]],2,sum))), nst,it, dimnames = list(stnms, 1:it))  # [nst,it]
                            
    TAC <- ifelse(B*CT < TAC.yr*QS.ss, B*CT, TAC.yr*QS.ss) 
    
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

    if(fleets.ctrl[[flnm]]$restriction == 'catch'){

        efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
        effs <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))
        Cr.f <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))

        q.m <- alpha.m <- beta.m  <- vector('list', length(sts))
        names(q.m) <- names(alpha.m) <- names(beta.m) <- sts 

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

            for(mt in mtnms){

                 if(!(st %in% names(fl@metiers[[mt]]@catches))) next
                    
                 q.m[[st]][mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE] 
                 alpha.m[[st]][mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE] 
                 beta.m[[st]][mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE] 
            }    
        
            Cr.f[st,] <- TAC[st,]*QS[[st]][flnm,]
                
            for(i in 1:it){          
                effort.fun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'effort', sep = '.')
                effs[st, i] <-  eval(call(effort.fun, Cr = Cr.f[st,i], B = B[st,i], Ba = Ba[[st]][,,,,,i,drop=F], q.m = q.m[[st]][,,,i,drop=F], 
                                efs.m = efs.m[,i], alpha.m = alpha.m[[st]][,,,i,drop=F], beta.m = beta.m[[st]][,,,i,drop=F]))
            }
        }
            
        # Choose the effort.
        eff <- effRule.SMFB(effs = effs, prev.eff = matrix(fl@effort[,yr-1,,ss,drop=T],1,it), 
                            rule = fleets.ctrl[[flnm]]$effort.restr)
            
        # Capacity restrictions.  
        eff <- capacityRest.SMFB(eff, c(fl@capacity[,yr,,ss,drop=T]))                                   
        fleets[[flnm]]@effort[,yr,,ss] <- eff 
          
    #    save(advice,alpha.m,B,beta.m,Cr.f,CT,eff,effs,efs.m,fleets.ctrl, 
    #         q.m,QS,QS.ss,TAC,TAC.yr, file = paste(flnm, file = '.RData', sep = ""))
   
        # Update the quota share of this step and the next one if the 
        # quota share does not coincide with the actual catch. (update next one only if s < ns).
        for(st in sts){
         #   browser()
            yr.share       <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [it]
            ss.share       <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, it))# [it,ns]
            quota.share.OR <- matrix(t(yr.share*ss.share), ns, it)
            # The catch.
            catchFun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'CatchFleet', sep = ".")
            catch <- eval(call(catchFun, Ba = Ba[[st]], B = B[st,], effort = eff, efs.m = efs.m, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]]))
            
            quota.share    <- updateQS.SMFB(QS = quota.share.OR, TAC = TAC.yr[st,], catch = catch, season = ss)        # [ns,it]
                              
            fleets.ctrl$seasonal.share[[st]][flnm,yr,,] <- t(t(quota.share)/apply(quota.share, 2,sum)) #[ns,it], doble 't' to perform correctly de division between matrix and vector.
        }
    }
    
    else{ # landings restriction.
        stop('Not yet implemented')
    } 
    
    return(list(fleets = fleets, fleets.ctrl = fleets.ctrl))
}

#-------------------------------------------------------------------------------
# SSFB (fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1)
#-------------------------------------------------------------------------------
SSFB <- function(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1){

    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )

    # If year/season/iter numerics => indicate position
    # else names => get positions.

    # 'year' dimension.
    yr <- year
    if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
    if(length(yr) == 0) stop('The year is outside object time range')

    # 'season' dimension.
    ss <- season
    if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
    if(length(ss) == 0) stop('The season is outside object season range')

    # Check fleets.ctrl elements.
    if(!(fleets.ctrl[[flnm]]$restriction %in% c('catch', 'landings')))
        stop("fleets.ctrl[[flnm]]$restriction must be equal to 'catch' or 'landings'")

    # Dimensions.
    nf  <- length(fleets);         flnms <- names(fleets)
    nst <- length(biols);          stnms <- names(biols)
    ns  <- dim(biols[[1]]@n)[4]
    it  <- dim(biols[[1]]@n)[6]

    # Data
    # - seasonal biomass
    B    <- matrix(t(sapply(stnms, function(x){   # biomass in the middle if age struc. of the season  [nst,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss, drop=T])
                                else return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss, drop=T])})) , nst,it, dimnames = list(stnms, 1:it))

    # - stock quota share by fleet
    QS   <- lapply(stnms, function(x){           # list of stocks, each stock [nf,it]
                            # Calculate QS by fleet for the year and season
                            yr.share    <- advice$quota.share[[x]][,yr,, drop=T]        # [nf,it]
                            ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss, drop=T]   # [nf,it]
                            quota.share <-  matrix(yr.share*ss.share, nf, it, dimnames = list(flnms, 1:it))
                            quota.share[is.na(quota.share)] <- 0
                            return(quota.share)})
    names(QS) <- stnms

    # If TAC >= B*alpha => TAC = B*alpha.
    # - annual TAC for each stock
    TAC.yr  <- matrix( advice$TAC[,yr,drop=T], nst, it, dimnames=list(stnms,1:it))  # [nst,it]
    # - limit percentage of stock to be captured
    CT      <- fleets.ctrl$catch.threshold[,yr,,ss, drop=T]  # [nst,it]
    # - stock quota share (sum of all fleets)
    QS.ss   <- matrix(t(sapply(stnms, function(x) apply(QS[[x]],2,sum))), nst,it, dimnames = list(stnms, 1:it))  # [nst,it]

    # - seasonal TAC per stock
    TAC <- ifelse(B*CT < TAC.yr*QS.ss, B*CT, TAC.yr*QS.ss)

    # Re-scale QS to fleet share within the season instead of season-fleet share within year.
    QS   <- lapply(stnms, function(x){          # list of stocks, each stock [nf,it]
                            res <- sweep(QS[[x]], 2, apply(QS[[x]],2, sum), "/")
                            res[is.na(res)] <- 0
                            return(res)})
    names(QS) <- stnms
    
    fl    <- fleets[[flnm]]
    sts   <- catchNames(fl)
    mtnms <- names(fl@metiers)

    if(fleets.ctrl[[flnm]]$restriction == 'catch'){

        st_mt <- character(length(sts)); names(st_mt) <- sts
        for(st in sts){
          # Search the metier targeting st, check that there is only one
          st.mt <- c()
          for (mt in mtnms) if( st %in% names(fl@metiers[[mt]]@catches) ) st.mt <- c(st.mt,mt) 
          if (length(st.mt)>1)
            stop( paste("There is a metier targeting more than one stock, therefore '",
                          flnm,"' fleet is not valid to perform SSFB",sep="")) 
          st_mt[st] <- st.mt
        }

        # Effort by metier (therefore by stock)
        effs.m <- matrix(NA, length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
        
        q.m <- alpha.m <- beta.m <- array(0, dim = c(length(mtnms), it, length(sts)), dimnames = list(metier = mtnms, 1:it, stock = sts))
        
        for(st in sts){     # q.m, alpha.m.... by metier (which always is stock specific)
          q.m[st_mt[st],,st]     <- fl@metiers[[st_mt[st]]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE] 
          alpha.m[st_mt[st],,st] <- fl@metiers[[st_mt[st]]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE] 
          beta.m[st_mt[st],,st]  <- fl@metiers[[st_mt[st]]]@catches[[st]]@beta[,yr,,ss, drop = TRUE]
          
          for(i in 1:it)
            effs.m[st_mt[st], i] <- cobbDougInv(Cr = TAC[st,i]*QS[[st]][flnm,i], B = B[st,i], q.m = q.m[st_mt[st],i,st],
                                          alpha.m = alpha.m[st_mt[st],i,st], beta.m = beta.m[st_mt[st],i,st])
        }

        # Capacity restrictions        
        # (i.e. if total fleet effort higher than effective capacity, then effort reduction)
        nmt <- length(st_mt)
        rule <- fleets.ctrl[[flnm]]$effReduce.rule
        effDay.perc <- fleets.ctrl[[flnm]]$effectiveDay.perc
        if (is.null(rule)) { # reduce the same proportion in all the metiers
            effs.m <- capacityRest.SSFB( eff.m=effs.m, capacity=c(fl@capacity[,yr,,ss,drop=T]*effDay.perc[,yr,,ss,drop=T]))
        } else if ( rule == 'month.price' ) { # we use the previous year price as weighting factor in effort reduction
            prices <- matrix(t(sapply(sts, function(x)
                                              apply(fl@metiers[[st_mt[x]]]@catches[[x]]@price[,yr-1,,ss,],c(2,4:6),mean)[drop=T])), 
                                length(sts), it, dimnames = list(sts, 1:it))
            weight.price <- prices; rownames(weight.price) <- st_mt[rownames(weight.price)]
            for (i in 1:it) for (j in 1:length(st_mt)) if (is.na(weight.price[j,i])) weight.price[j,i] <- 0 
            price.order <- matrix( , length(st_mt), it, dimnames=list(c(), 1:it))
            for (i in 1:it) price.order[,i] <- names(sort(weight.price[,i],decreasing=T)) 
            effs.m <- capacityRest.SSFB( eff.m=effs.m, capacity=c(fl@capacity[,yr,,ss,drop=T]*effDay.perc[,yr,,ss,drop=T]), order.m=price.order)
        } else
            stop(paste("Invalid rule to reduce metier efforts in fleet '",flnm,"'",sep=""))
        
        # Total fleet effort:
        fleets[[flnm]]@effort[,yr,,ss,] <- apply(effs.m,2,sum)
        # Effort proportion by metier:
        for (mt in st_mt)
          fleets[[flnm]]@metiers[[mt]]@effshare[,yr,,ss,] <- effs.m[mt,] / fleets[[flnm]]@effort[,yr,,ss,]
          
        eff <- fleets[[flnm]]@effort[,yr,,ss,drop=T] 
        efs.m <- matrix(t(sapply(mtnms, function(x) fleets[[flnm]]@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))

        # Update the quota share of this step and the next one if the 
        # quota share does not coincide with the actual expected catch. (update next one only if s < ns).
        for(st in sts){
            yr.share       <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [it]
            ss.share       <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, it))# [it,ns]
            quota.share.OR <- matrix(t(yr.share*ss.share), ns, it)
            quota.share    <- updateQS.SSFB(QS = quota.share.OR, TAC = TAC.yr[st,], B = B[st,], E = eff,          # [ns,it]
                                efs.m = efs.m, q.m = q.m[,,st], alpha.m = alpha.m[,,st], beta.m = beta.m[,,st], season = ss)
            fleets.ctrl$seasonal.share[[st]][flnm,yr,,] <- t(t(quota.share)/apply(quota.share, 2,sum)) #[ns,it], doble 't' to perform correctly de division between matrix and vector.
        }

    } else{ # landings restriction.
        stop('Not yet implemented')
    }
    
    return(list(fleets = fleets, fleets.ctrl = fleets.ctrl))
}
