#-------------------------------------------------------------------------------
#       Update landings, discards, landings.n, discards.n slots.
#
# - aggregated.CobbDoug: Catch at age calculation when Cobb douglas is applied at
#       total biomass level.
# - ageBased.CobbDoug: Catch at age calculation when Cobb douglas is applied at
#       age level.
#  - updateCatch: Update the slots related to catch using the appropiate function.
# 
# Dorleta GarcYYYa
# Created: 28/10/2010 12:33:04
# Changed:03/06/2011 07:53:53
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# updateCatch(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
updateCatch <- function(fleets, biols, advice, fleets.ctrl, year = 1, season = 1){

    fleet.names <- names(fleets)
    
    fleets <- unclass(fleets) # convert it into a list to speed up de computations
    
    for(flnm in fleet.names){
        # Which stocks are caught by fleet flnm.
        flsts <- catchNames(fleets[[flnm]])
        for(st in flsts){
            catch.model <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'CAA', sep = ".")
            fleets <- eval(call(catch.model, fleets = fleets, biols = biols, fleets.ctrl = fleets.ctrl, advice = advice, year = year, season = season, flnm = flnm, stknm = st))
        }
    }
    
     fleets <- FLFleetsExt(fleets)
    
    # Correct the catch in case Ca > Ba Or C > B
    fleets <- CorrectCatch(fleets = fleets, biols = biols, year = year, season = season)
    
    
    return(fleets)
}

#-------------------------------------------------------------------------------
# CobbDouglasBio.CatchFleet(effort, Ba, q.m, efs.m, alpha.m, beta.m)
#-------------------------------------------------------------------------------
CobbDouglasBio.CatchFleet <- function(effort, Ba, B, q.m, efs.m, alpha.m, beta.m,...){

    nmt <- dim(efs.m)[1]
    it  <- dim(efs.m)[2]
    
   ## Redimensionate all the objects into dimension [nmt,it]
    
    # dim(q.m) = dim(alpha.m) = dim(beta.m) = [nmt,na,nu,it]
    q.m     <- matrix(q.m[,,,,drop=TRUE],nmt,it)      # [nmt,it]
    alpha.m <- matrix(alpha.m[,,,,drop=TRUE],nmt,it)  # [nmt,it]
    beta.m  <- matrix(beta.m[,,,,drop=TRUE],nmt,it)   # [nmt,it]

    # dim(B) = dim(effort) = [it]
    B       <- matrix(B, nmt, it, byrow = TRUE)      # [nmt,it]
    effort  <- matrix(effort, nmt, it, byrow = TRUE) # [nmt,it]

    catch <- apply(q.m*(B^beta.m)*((effort*efs.m)^alpha.m),2,sum) # sum catch along metiers
    
    return(catch)  
}


#-------------------------------------------------------------------------------
# CobbDouglasAge.CatchFleet(effort, Ba, q.m, efs.m, alpha.m, beta.m)
#-------------------------------------------------------------------------------
CobbDouglasAge.CatchFleet <- function(effort, Ba, q.m, efs.m, alpha.m, beta.m,...){

    Ba      <- Ba[,,,,,,drop=T]         # [naYYY?,nuYYY?,itYYY?]   nu and(or it can be equal 1 => the dimension would be dropped
    
    dimq  <- dim(q.m)
    zz    <- ifelse(dimq == 1, FALSE, TRUE)
    
    Ba <- array(Ba, dim = c(dim(q.m)[2:4], dim(q.m)[1]))  # [na,nuYYY?,itYYY?,mt]
    Ba <- aperm(Ba, c(4,1:3))  # [mt,na,nu,it]
    
    effort <- array(effort, dim = c(length(effort), dim(q.m)[1:3])) # [it,mt,na,nu]
    effort <- aperm(effort, c(2:4,1))  # [mt,na,nu,it]
    
    efs.m <- array(efs.m, dim = c(dim(efs.m), dimq[2:3]))
    efs.m <- aperm(efs.m, c(1,3:4,2))
    
    catch <- apply(q.m*(Ba^beta.m)*((effort*efs.m)^alpha.m), 4,sum)
    
    return(catch)  
}



#-------------------------------------------------------------------------------
# aggregated.CobbDoug(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
CobbDouglasBio.CAA  <- function(fleets, biols, fleets.ctrl, advice, year = 1, season = 1, flnm = 1, stknm = 1, ...){

    nf    <- length(fleets)
    stnms <- names(biols)
    nst   <- length(stnms)
    it    <- dim(biols[[1]]@n)[6]
    
    yr <- year
    ss <- season
    f  <- flnm
    st <- stknm
    
    fleets <- unclass(fleets)

    fl    <- fleets[[f]]
    sts   <- catchNames(fl)
    mtnms <- names(fl@metiers)
    
    if(!(st %in% sts)) return(fleets)
                                             
    tac <- rep('Inf',it)
    
    # if TAC overshoot is discarded, calculate seasonal TAC to calculate the discards.
    TACOS <- fleets.ctrl[[flnm]][[stknm]][['discard.TAC.OS']]
    TACOS <- ifelse(is.null(TACOS), TRUE, TACOS) 
    if(TACOS){
        yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
        ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
        QS          <- yr.share*ss.share                                          # [it]
        QS[is.na(QS)] <- 0              
        tac <- (advice$TAC[st,yr]*QS)[drop=T] # it
    }
    
    # biomass in the middle if age struc. of the season  B[it]
    B <- ifelse(dim(biols[[st]]@n)[1] > 1,
                    unitSums(quantSums(biols[[st]]@n*biols[[st]]@wt*exp(-biols[[st]]@m/2)))[,yr,,ss, drop=T],
                    (biols[[st]]@n*biols[[st]]@wt)[,yr,,ss, drop=T])
    
    Ba <- biols[[stknm]]@n[,yr,,ss]*biols[[stknm]]@wt[,yr,,ss]*exp(-biols[[stknm]]@m[,yr,,ss]/2)  # Ba[na,1,1,1,1,it]
            
    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    eff   <- matrix(fl@effort[,yr,,ss],1, it, dimnames = list(1, 1:it))
                     
    
    for(mt in mtnms){
       #       print(c(f, " - ", st, " - ", mt))
        if(!(st %in% names(fl@metiers[[mt]]@catches))) next

            cobj <- fl@metiers[[mt]]@catches[[st]]

            q.m         <- cobj@catch.q[,yr,,ss, drop = TRUE]
            alpha.m     <- cobj@alpha[,yr,,ss, drop = TRUE]
            beta.m      <- cobj@beta[,yr,,ss, drop = TRUE]

            Ctotal <- array((eff*efs.m[mt,])^alpha.m*B^beta.m*q.m,dim = c(rep(1,5),it))

            tac.disc <- ifelse(Ctotal[,,,,,,drop=T] < tac, 1, tac/Ctotal[,,,,,,drop=T])

            dsa <- cobj@discards.sel[,yr,,ss]  # [na,1,nu,1,1,it]
            lsa <- cobj@landings.sel[,yr,,ss]  # [na,1,nu,1,1,it]
            sa  <- (dsa + lsa)  
            
            # Recalculate dsa and lsa according to 'tac.disc'     # [na,nu,it]
            lsa <- lsa*tac.disc
            dsa <- sa - lsa

            if(dim(biols[[st]]@n)[1] == 1){
                cobj@discards[,yr,,ss]   <- Ctotal*dsa/(sa*tac.disc)
                cobj@landings[,yr,,ss]   <- Ctotal*lsa*tac.disc/sa
                cobj@discards.n[,yr,,ss] <- cobj@discards[,yr,,ss]/cobj@discards.wt[,yr,,ss]
                cobj@landings.n[,yr,,ss] <- cobj@landings[,yr,,ss]/cobj@landings.wt[,yr,,ss]
          
                fl@metiers[[mt]]@catches[[st]] <- cobj
            }
            else{ # age structured stock, biomass level Cobb Douglas.
            #     browser() 
                Bs <- apply(sa*Ba, 6,sum)           # it
                Ca <- sweep(sa*Ba, c(2,4:6), Ctotal/Bs, "*") # [na,nu,it]

                cobj@discards.n[,yr,,ss] <- Ca*dsa/sa/cobj@discards.wt[,yr,,ss]
                cobj@landings.n[,yr,,ss] <- Ca*lsa/sa/cobj@landings.wt[,yr,,ss]

                # When sa = 0 <-  land.n & dis.n = NA => change to 0.
                cobj@landings.n[,yr,,ss][is.na(cobj@landings.n[,yr,,ss])] <- 0
                cobj@discards.n[,yr,,ss][is.na(cobj@discards.n[,yr,,ss])] <- 0

                cobj@discards[,yr,,ss] <- apply(Ca*dsa/sa,c(2,4,6),sum)
                cobj@landings[,yr,,ss] <- apply(Ca*lsa/sa,c(2,4,6),sum)

                fl@metiers[[mt]]@catches[[st]] <- cobj
            }
        }
          
    fleets[[f]] <- fl
    
    fleets <- FLFleetsExt(fleets)
      
    return(fleets)
}

            
#-------------------------------------------------------------------------------
# ageBased.CobbDoug(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
CobbDouglasAge.CAA <- function(fleets, biols, fleets.ctrl, advice, year = 1, season = 1, flnm = 1, stknm = 1,...){

    nf    <- length(fleets)
    stnms <- names(biols)
    nst   <- length(stnms)
    it    <- dim(biols[[1]]@n)[6]
    
    fleets <- unclass(fleets)
    
    yr <- year
    ss <- season
    st <- stknm
    
    fl    <- fleets[[flnm]]
    sts   <- catchNames(fl)
    mtnms <- names(fl@metiers)

    if(!(st %in% sts)) return(fleets)
    
    tac <- rep('Inf',it)
    
    # if TAC overshoot is discarded, calculate seasonal TAC to calculate the discards.
    TACOS <- fleets.ctrl[[flnm]][[stknm]][['discard.TAC.OS']]
    TACOS <- ifelse(is.null(TACOS), TRUE, TACOS) 
    if(TACOS){
        yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
        ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
        QS          <- yr.share*ss.share                                          # [it]
        QS[is.na(QS)] <- 0              
        tac <- (advice$TAC[st,yr]*QS)[drop=T] # it
    }
    
    if(dim(biols[[st]]@n)[1] == 1) stop(st, ' stock has no ages, Cobb Douglas cannot be applied at age level then! correct the "catch.model" argument in "fleets.ctrl" argument!\n')
    
    Ba <- (biols[[st]]@n*biols[[st]]@wt*exp(-biols[[st]]@m/2))[,yr,,ss]  # Ba[na,it], biomass at age in the middle  of the season,

    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                        length(mtnms), it, dimnames = list(metier = mtnms, 1:it))    # [nmt,it]
    eff   <- matrix(fl@effort[,yr,,ss],1, it, dimnames = list(1, 1:it))              # [1,it]
                 
                 
    for(mt in mtnms){
        
     #       print(c(f, " - ", st, " - ", mt))
        if(!(st %in% names(fl@metiers[[mt]]@catches))) next
            
        cobj <- fl[[mt]][[st]]

        q.m         <- cobj@catch.q[,yr,,ss]     # [na,1,nu,1,1,it] : na & nu  can be 1.
        alpha.m     <- cobj@alpha[,yr,,ss]       # [na,1,nu,1,1,it] : na & nu  can be 1.
        beta.m      <- cobj@beta[,yr,,ss]        # [na,1,nu,1,1,it] : na & nu  can be 1.

        na <- dim(q.m)[1]
        nu <- dim(q.m)[3]

        efm <- array(eff*efs.m[mt,], dim = c(it,na,1,nu,1,1))
        efm <- aperm(efm, c(2:6,1))

        Ca <- q.m*efm^alpha.m*Ba^beta.m   # [na,1,nu,1,1,it] : na & nu  can be 1.

        Ct <- apply(Ca, 6, sum)
        
        tac.disc <- ifelse(Ct < tac, 1, tac/Ct)

        dsa <- cobj@discards.sel[,yr,,ss]  # [na,1,nu,1,1,it]
        lsa <- cobj@landings.sel[,yr,,ss]  # [na,1,nu,1,1,it]
        sa  <- (dsa + lsa)  
            
        # Recalculate dsa and lsa according to 'tac.disc'     # [na,nu,it]
        lsa <- lsa*tac.disc
        dsa <- sa - lsa             # [na,nu,it]

        cobj@discards.n[,yr,,ss] <- Ca*dsa/sa/cobj@discards.wt[,yr,,ss]
        cobj@landings.n[,yr,,ss] <- Ca*lsa/sa/cobj@landings.wt[,yr,,ss]

        # When sa = 0 <-  land.n & dis.n = NA => change to 0.
        cobj@landings.n[,yr,,ss][is.na(cobj@landings.n[,yr,,ss])] <- 0
        cobj@discards.n[,yr,,ss][is.na(cobj@discards.n[,yr,,ss])] <- 0

        cobj@discards[,yr,,ss] <- apply(Ca*dsa/sa,c(2,4,6),sum,na.rm=T)
        cobj@landings[,yr,,ss] <- apply(Ca*lsa/sa,c(2,4,6),sum,na.rm=T)

        fl@metiers[[mt]]@catches[[st]] <- cobj          
    }
    
    fleets[[flnm]] <- fl
    
    fleets <- FLFleetsExt(fleets)
    
    return(fleets)
}
            
            
#-------------------------------------------------------------------------------
# CorrectCatch(fleets, biols, year = 1, season = 1)
# Given that in some production functions it can happen that
# Ca > Ba (age struc. pop) or C > B (bio struc. pop), if this happens the
# catch is corrected setting Ca = Ba and C = B. For this end the catch is reduced
# in the same degree in all the fleets. This could imply a 'revision' in the
# production function parameters for the year, season and iteration in question
# but it is not done internally because it does not affect other steps and it
# would imply a los in the generality of the 'CorrectCatch' function because in
# for doing so the catch production function should be used.
#-------------------------------------------------------------------------------

CorrectCatch <- function(fleets, biols, year = 1, season = 1,...){

  #  fleets <- unclass(fleets)
    yr <- year
    ss <- season
    it    <- dim(biols[[1]]@n)[6]
    nst <- length(biols)
    
    stnms <- names(biols)
    flnms <- names(fleets)

    Ba   <- lapply(stnms, function(x){   # biomass at age in the middle  of the season, list elements: [na,it] if age structured, [1,it] if biomass.
                            if(dim(biols[[x]]@n)[1] > 1)
                                return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss])
                            else return(matrix((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss],1,it))})
    names(Ba) <- stnms

    B    <- matrix(t(sapply(stnms, function(x){   # biomass in the middle if age struc. of the season  [ns,it]
                if(dim(biols[[x]]@n)[1] > 1)
                    return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss, drop=T])
                else return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss, drop=T])})) , nst,it, dimnames = list(stnms, 1:it))


    for(st in stnms){
   #    print(st)
#   if(st == 'CANK') browser()
        if(dim(Ba[[st]])[1] > 1){ # age structured

            Cat  <- catchStock(fleets, st)[,yr,,ss]
            
            # Convert the [age,unit] combination into a continuous age.
            Ba.  <- unit2age(Ba[[st]]) # [na*nu,1,1,it]
            Cat. <- unit2age(Cat)
            K.   <- array(1,dim = dim(Ba.))   # Catch multipliers

            # CORRECT Ca if Ca > Ba, the correction is common for all the fleets.
            for(i in 1:it){

                if(any((Ba[[st]][,,,,,i] - Cat[,,,,,i]) < 0)){

                    cat('Ba < Ca, for some "a" in stock',st, ', and iteration ', i,  '\n')

                    a.minus         <- which(Ba.[,,,i] < Cat.[,,,i])
                    a.plus          <- which(Ba.[,,,i] >= Cat.[,,,i])
                    K.[a.minus,,,i] <- 0.99*Ba.[a.minus,,,i]/Cat.[a.minus,,,i]

                    # The correction below would correspond with a compensation of the decrease in 'a.minus' ages.
                    # K.[a.plus,,,i]  <- (Ct[i] - sum(Ba.[a.minus,,,i]))/sum(Cat.[a.plus,,,i])
                }
            }
            K <- age2unit(K., Ba[[st]])  # [na,1,nu,1,1,it]
            K[1,,-(1:ss)] <- 1   # This recruits do not exists  yet.
        }
        else{ # biomass dynamic.
        #   browser()
            Ct  <- catchWStock(fleets, st)[,yr,,ss,drop = F]
            Bst <- array(B[st,], dim = c(1,1,1,1,1,it))     # [it]
            K   <- rep(1,it)  # Catch multipliers

            if(any((Bst[,,,,,i] - Ct[,,,,,i]) < 0)){
                i.minus <-  which((Bst[,,,,,i] - Ct[,,,,,i]) < 0)
                cat('B < C, for  stock',st, ', and iteration(s) ', i.minus,  '\n')
                K[i.minus] <- Bst[i.minus]/Ct[i.minus]

            }

            K <- FLQuant(K, dimnames = dimnames(biols[[st]]@n[,yr,,ss]))
        }

        if(all(K==1)) next
        # Correct the catch fleet by fleet.
        # does the fleet_metier catch the stock??  which of them catch the stock??
        flinfo  <- stock.fleetInfo(fleets)
        flmtpos <-  colnames(flinfo)[which(flinfo[st, ] == 1)]

        for(k in flmtpos){
            k. <- strsplit(k, "&&")[[1]]
            fl <- k.[1]
            mt <- k.[2]

            cobj <- fleets[[fl]][[mt]][[st]]
            
     #       print(c(cobj@landings.n[14,yr,,ss]))

            cobj@landings.n[,yr,,ss] <-  cobj@landings.n[,yr,,ss]*K
            cobj@discards.n[,yr,,ss] <-  cobj@discards.n[,yr,,ss]*K
            cobj@landings[,yr,,ss]   <-  apply((cobj@landings.n*cobj@landings.wt)[,yr,,ss], c(2,4,6), sum)
            cobj@discards[,yr,,ss]   <-  apply((cobj@discards.n*cobj@discards.wt)[,yr,,ss], c(2,4,6), sum)

     #        print(c(cobj@landings.n[14,yr,,ss]))
             
           fleets[[fl]]@metiers[[mt]]@catches[[st]] <- cobj
        }
    }
    
#    fleets <- FLFleetsExt(fleets)
    return(fleets)
}




                            
                
        
            
             

    
        
    

