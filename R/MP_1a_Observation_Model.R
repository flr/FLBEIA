#-------------------------------------------------------------------------------
#                       OBSERVATION MODEL FUNCTIONS
#
#   - perfectObs: Variables in FLStock are observed without error.  
#   - age2ageDat: Create an age structured FLStock from an age structured FLBiol. 
#              (Aging error, underreporting....).Only the data neccesary for the 
#              assessment is generated (NO stock.n, NO harvest)
#              slightly modified from rroa's original functions.
#   - age2agePop: Create an age structured FLStock from an age structured FLBiol. 
#              (Aging error, underreporting....). stock.n and harvest _are_ observed
#              slightly modified from rroa's original functions.
#   - bio2bioDat: Create an aggregated FLStock from an aggregated FLBiol. 
#              (underreporting....). Only the data neccesary for the 
#              assessment is generated (NO stock, NO harvest)
#              slightly modified from rroa's original functions.
#   - bio2bioPop: Create an aggregated FLStock from an aggregated FLBiol. 
#              (underreporting....). As bio2bioDat but stock and harvest
#               _are_ observed.
#              slightly modified from rroa's original functions.
#   - age2bioDat: Create an aggregated FLStock from an age structured FLBiol. 
#              (underreporting....). Only the data neccesary for the 
#              assessment is generated (NO stock, NO harvest)
#              slightly modified from rroa's original functions.
#   - age2bioPop: Create an aggregated FLStock from an age structured FLBiol. 
#              (underreporting....). As age2bioDat but stock and harvest
#               _are_ observed.
#              slightly modified from rroa's original functions.
#   - ageInd: Update an age structured FLIndex from an age structured FLBiol 
#             (aging error + multiplicative error)  (dga)
#   - bioInd: Update a FLIndex aggregated in biomass from an age structured or 
#             biomass aggregated FLBiol (multiplicative error)  (dga)
#   - ssbInd: Update a FLIndex aggregated in biomass from an age structured FLBiol and a FLFleet object
#             (ssb + multiplicative error)  (ssm)
#   - cbbmInd: Update a FLIndex aggregated in biomass (age 1, 2+) from an age structured FLBiol and a FLFleets object
#             (B1 and B2+ + multiplicative error)  (ssm)
#
# Dorleta GarcYYYa
# Created: 09/12/2010 14:37:43
# Changed: 15/04/2013 10:59:56 (correct a bug in PerfectObs problem with seasons) 
#(rroa's functions inserted)
#         2013-06-07 12:20:04 Sonia Sanchez - code revised and names changed for coherence
#         2014-02-24 17:47:21 Sonia Sanchez - ssbInd and cbbmInd functions added
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# perfectObs(biol, fleets, covars, obs.ctrl, year = 1, season = 1)
#-------------------------------------------------------------------------------
perfectObs <- function(biol, fleets, covars, obs.ctrl, year = 1, season = NULL, ...){

    # THE ASSESSMENT IS BEING CARRIED OUT IN <year> => THE 'OBSERVATION' GOES UP TO <year-1>
    
    st <- biol@name
    na <- dim(biol@n)[1]
    ns <- dim(biol@n)[4]
    it <- dim(biol@n)[6]
    ss <- ifelse(is.null(season), dim(biol@n)[4], season)
    
    if ( year > dims(biol)$year) biol <- window( biol, start=dims(biol)$minyear, end=dims(biol)$maxyear+1)
    
    # FIRST SEASON, FIRST UNIT:
    # biol@wt = "stock.wt" = "catch.wt" = "discards.wt" = "landings.wt" = "mat" = "harvest.spwn" = "m.spwn
    res <- propagate(as(biol, 'FLStock')[,1:(year-1),1,1], it, fill.iter = TRUE)
    dimnames(res) <- list(unit="unique")
    
    #res@range[c(1:3,6:7)] <- biol@range[c(1:3,6:7)]
    #names(res@range[6:7]) <- c('minfbar', 'maxfbar')
        
    res@discards.wt <- res@landings.wt <- res@catch.wt <- res@stock.wt
        
    # "stock.n":  FIRST SEASON and SUM ALONG UNITS except recruitment
    # rec = n[1,,1,1] + n[1,,2,2] + n[1,,3,3] + n[1,,4,4]
    # n up to (year) to use it after in the 'f' calculation.
    n <- unitSums(biol@n)[,1:year,,1]
    n[n == 0] <- 1e-6   # if n == 0 replace it by a small number to avoid 'Inf' in harvest.
        
    if(dim(biol@n)[3] > 1){
        for(u in 2:dim(biol@n)[3])
            n[1,] <- n[1,] + biol@n[1,1:year,u,u]
    } else {
        for (s in c(1:ns)) { 
          n[1,1:(year-1),] <- biol@n[1,1:(year-1),1,s,]
          if( sum( n[1,1:(year-1),] != 0, na.rm = T) > 0) break # spawning season
        }  
    }
    # for current year if season before recruitment season:
    if (ss != ns)
      n[1,1:(year-1),] <- ifelse( is.na(n[1,1:(year-1),]), 0, n[1,1:(year-1),])
         
    stock.n(res) <- n[,1:(year-1)]
    
    stock(res) <- quantSums(res@stock.n*res@stock.wt)
        
    # SUM ALONG SEASONS AND FIRST UNIT: "m"
    m(res)[]      <- seasonSums(biol@m)[,1:(year-1),1,]
    m.spwn(res)[] <- seasonSums(spwn(biol))[,1:(year-1),1,]/ns
    if (ss < ns){ # sum only along s<=ss for last year
      m(res)[,year-1,]      <- seasonSums(biol@m[,year-1,1,1:ss,])
      m.spwn(res)[,year-1,] <- seasonSums(spwn(biol)[,year-1,1,1:ss])/length(1:ss)
    }

    # SUM ALONG UNITS AND SEASONS, OBTAINED FROM FLFLEETS: 
    # "catch", "catch.n", "discards"     "discards.n" "landings"     "landings.n"
    land.n <- apply(landStock(fleets, st), c(1:2,6),sum)[,1:(year-1),]
    disc.n <- apply(discStock(fleets, st), c(1:2,6),sum)[,1:(year-1),]
    dimnames(land.n)[1:5] <- dimnames(disc.n)[1:5] <- dimnames(landings.n(res))[1:5]
    landings.n(res) <- land.n
    discards.n(res) <- disc.n
    catch.n(res)    <- res@discards.n + res@landings.n
    landings(res)   <- quantSums(res@landings.n*res@landings.wt)
    discards(res)   <- quantSums(res@discards.n*res@discards.wt)
    catch(res)      <- quantSums(res@catch.n*res@catch.wt)
        
    # harvest: * if age structured calculate it from 'n'.
    #          * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
    if(na == 1){
        harvest(res)[] <- (res@catch)/(res@stock.n*res@stock.wt) 
        units(res@harvest) <- 'hr'
    } else{
        harvest(res) <- catch.n(res)*NA #! Artefact to avoid crashing (sets correct dim for iters)
        units(res@harvest) <- 'f'
        ai <- ifelse( ss==ns, na-1, na-2)
        res@harvest[-c(ai:na),] <- log(n[-c(ai:na),-year]/n[-c(1,(ai+1):na),-1]) - res@m[-c(ai:na),]
        # F < 0 -> correct F=0
        res@harvest[-c(ai,na),] <- ifelse( is.na(res@harvest[-c(ai,na),]) | as.numeric(res@harvest[-c(ai,na),])<0, 0, res@harvest[-c(ai,na),])
        
        n. <- array(res@stock.n[drop=T], dim = c(na,year-1,it))     # [na,ny,it]
        m. <- array(res@m[drop=T], dim = c(na,year-1,it))            # [na,ny,it]
        c. <- array(res@catch.n[drop=T], dim = c(na,year-1,it))      # [na,ny,it]
        
        
        fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
        
        for(y in 1:(year-1)){
          
            for(a in ai:na){
                for(i in 1:it){
                
                if(is.na(n.[a,y,i])) res@harvest[a,y,,,,i]<-0
                
                if(!is.na(n.[a,y,i])) {
                    n.[a,y,i] - c.[a,y,i]
                    if(n.[a,y,i] < c.[a,y,i]) res@harvest[a,y,,,,i] <- 10
                    
                    else{
                        zz <- try(ifelse(n.[a,y,i] == 0 | c.[a,y,i] == 0, 0,
                                                uniroot(fobj, lower = 1e-300, upper = 1e6, n = n.[a,y,i], m=m.[a,y,i], c = c.[a,y,i])$root), silent = TRUE)  
                        res@harvest[a,y,,,,i] <- ifelse(is.numeric(zz), zz, res@harvest[ai-1,y,,,,i] )
                    }
                }
                }
            }
        }
        # for current year if season before recruitment season:
        if (ss != ns)
          res@harvest[1,year-1,] <- ifelse( is.na(res@harvest[1,year-1,]), 0, res@harvest[1,year-1,])
        
        ctot.age <- apply(landStock(fleets, st), c(1:2,4,6),sum)[,1:(year-1),] + apply(discStock(fleets, st), c(1:2,4,6),sum)[,1:(year-1),]
        ctot     <- seasonSums(ctot.age)
        c.perc <- ctot.age * NA
        for (s in c(1:ns)) c.perc[, , , s, ] <- ifelse(ctot==0,0, ctot.age[, , , s, ]/ctot)
        biol.spwn <- unitMeans(spwn(biol)[,1:(year-1),,])
        harvest.spwn(res)[] <- seasonSums(c.perc[,1:(year-1),,]*biol.spwn)

    }
    
    # If catc
    return(res)
}
 

#-------------------------------------------------------------------------------
# age2age(biol, fleets, advice, obs.ctrl, year, stknm)
# Age-Structured Observation of age structured pop
#  ** obs.ctrl in this case is a subset of the original obs.ctrl
#       obs.ctrl <- obs.ctrl[[stknm]][['stkObs']] when calling to age2age in obs.model function.
#-------------------------------------------------------------------------------
age2ageDat <- function(biol, fleets, advice, obs.ctrl, year, stknm,...){
    
    yr                <- year
                                 
    na                <- dim(biol@n)[1]
    ny                <- yr - 1
    it                <- dim(biol@n)[6]
    
    if (is.null(obs.ctrl$TAC.ovrsht))
      stop(paste("'TAC.ovrsht' array not defined for stock '",stknm,"'",sep=""))
    # TAC.ovrsht  can be numeric with dimension [1]  or FLQuant with dimension [1,dim(biol@n)[2],1,1,1,it]     
    # If TAC.ovrsht is numeric => convert it into an FLQuant. 
    if(is.null(dim(obs.ctrl$TAC.ovrsht))) obs.ctrl$TAC.ovrsht <- FLQuant(obs.ctrl$TAC.ovrsht, dim = c(1,dim(biol@n)[2],1,1,1,it))
                                                                
    ages.error        <- obs.ctrl$ages.error
    nmort.error       <- obs.ctrl$nmort.error[,1:ny,,drop=F]
    fec.error         <- obs.ctrl$fec.error[,1:ny,,drop=F]
    land.wgt.error    <- obs.ctrl$land.wgt.error[,1:ny,,drop=F]
    disc.wgt.error    <- obs.ctrl$disc.wgt.error[,1:ny,,drop=F]
    land.nage.error   <- obs.ctrl$land.nage.error[,1:ny,,drop=F]
    disc.nage.error   <- obs.ctrl$disc.nage.error[,1:ny,,drop=F]   
    TAC.ovrsht        <- obs.ctrl$TAC.ovrsht[,1:ny,,drop=F]
        
    if(is.null(ages.error)){
        ages.error <- array(0,dim = c(na, na, ny,it))
        for(a in 1:na) ages.error[a,a,,] <- 1
    }
    
    if(dim(ages.error)[1] != na | dim(ages.error)[2] != na)
         stop("ages.error array must have dim[1:2] identical to number of ages in stock")
     if(any(round(apply(ages.error,c(1,3:4), sum),2) != 1))
         stop("Some rows in ages.error array  don't add up to 1")
    
    for (e in c('nmort.error', 'land.wgt.error', 'disc.wgt.error', 
                 'fec.error', 'land.nage.error', 'disc.nage.error')) {
      err <- get(e)
      if (is.null(err))
        stop(paste("'",e,"' array not defined for stock '",stknm,"'",sep=""))
      if(dim(err)[1] != na)
        stop(paste("'",e,"' array, for stock '",stknm,"', must have dim[1] identical to number of ages in stock",sep=""))
      if (sum(err<=0)>0 | sum(is.na(err))>0)
        stop(paste("check values in '",e,"' array for stock '",stknm,"' (required values > 0)",sep=""))
      
    }
         
    stck              <- propagate(as(biol, "FLStock")[,1:ny,1,1],it, fill.iter = TRUE)  
    dimnames(stck) <- list(unit="unique")

    landings.wt(stck)[] <- Obs.land.wgt(fleets, ages.error, land.wgt.error, yr, stknm)
    landings.n(stck)[]  <- Obs.land.nage(fleets, ages.error, land.nage.error, stck@landings.wt, yr, stknm)
    landings(stck)      <- quantSums(seasonSums(stck@landings.n*stck@landings.wt))
    
    # In landings.wt the error due to age depends on landings.n, but that on m and mat only depends on error itself
    # because it is suppose that the biological sampling is independent.
    
    m(stck)[]          <- Obs.nmort(biol, ages.error, nmort.error, yr)
    mat(stck)[]        <- Obs.mat(biol, ages.error, fec.error, yr)

    # compare the landings with the advice and depending on TAC.ovrsht report landings. the misresporting is 
    # reported homogeneously.
    landings(stck)[]  <- FLQuant(ifelse(unclass(stck@landings) > TAC.ovrsht*advice$TAC[stknm,1:ny], TAC.ovrsht*advice$TAC[stknm,1:ny], stck@landings))
    
    ovrsht.red        <- stck@landings/quantSums(unitSums(seasonSums(stck@landings.n*stck@landings.wt)))  # [1,ny,,,,it]
    
    ovrsht.red[is.na(ovrsht.red)] <- 1
     
    landings.n(stck)[]   <- sweep(stck@landings.n, 2:6, ovrsht.red, "*") #distributing the overshoot subreporting of bulk landings in biomass equally over ages
   
    discards.wt(stck)[]  <- Obs.disc.wgt(fleets, ages.error, disc.wgt.error,  yr, stknm)
    stck@discards.wt[is.na(stck@discards.wt)] <- stck@landings.wt[is.na(stck@discards.wt)]
    discards.n(stck)[]   <- Obs.disc.nage(fleets, ages.error, disc.nage.error, stck@discards.wt, yr, stknm)
    discards(stck)       <- quantSums(seasonSums(stck@discards.n*stck@discards.wt))
    
    catch(stck)        <- stck@landings + stck@discards
    catch.n(stck)      <- stck@landings.n + stck@discards.n
    catch.wt(stck)     <- (stck@landings.n*stck@landings.wt + stck@discards.n*stck@discards.wt)/(stck@landings.n + stck@discards.n)
    
 #   stck@harvest      <- FLQuant(NA,dim=c(na,ny,1,1,1,it), dimnames=list(age=biol@range[1]:biol@range[2], year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))

    stck@harvest.spwn[] <- 0 # FLQuant(NA,dim=c(na,ny,1,1,1,it),dimnames=list(age=biol@range[1]:biol@range[2], year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))
    stck@m.spwn[]       <- 0 # FLQuant(NA,dim=c(na,ny,1,1,1,it),dimnames=list(age=biol@range[1]:biol@range[2], year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))

    return(stck)
}


#-------------------------------------------------------------------------------
# age2age(biol, fleets, obs.ctrl, year, stknm)
# Age-Structured Observation of age structured pop
#  ** obs.ctrl in this case is a subset of the original obs.ctrl
#       obs.ctrl <- obs.ctrl[[stknm]][['stkObs']] when calling to age2age in obs.model function.
#-------------------------------------------------------------------------------
age2agePop <- function(biol, fleets, advice, obs.ctrl, year, stknm,...){
                         
    ages.error        <- obs.ctrl$ages.error
    stk.nage.error    <- obs.ctrl$stk.nage.error
    stk.wgt.error     <- obs.ctrl$stk.wgt.error
    
    yr <- year 
                                 
    na                <- dim(biol@n)[1]
    ny                <- yr-1
    it                <- dim(biol@n)[6]
    
    if(is.null(ages.error)){
        ages.error <- array(0,dim = c(na, na, ny,it))
        for(a in 1:na) ages.error[a,a,,] <- 1
    }
    
    if(dim(ages.error)[1] != na | dim(ages.error)[2] != na)
         stop("ages.error array must have dim[1:2] identical to number of ages in stock")
    if(any(round(apply(ages.error,c(1,3:4), sum),2) != 1))
         stop("Some rows in ages.error array  don't add up to 1")
    
    for (e in c('stk.nage.error', 'stk.wgt.error')) {
      err <- get(e)
      if (is.null(err))
        stop(paste("'",e,"' array not defined for stock '",stknm,"'",sep=""))
      if(dim(err)[1] != na)
        stop(paste("'",e,"' array, for stock '",stknm,"', must have dim[1] identical to number of ages in stock",sep=""))
      if (sum(err<=0)>0 | sum(is.na(err))>0)
        stop(paste("check values in '",e,"' array for stock '",stknm,"' (required values > 0)",sep=""))
    }
         
    stck              <- age2ageDat(biol, fleets, advice, obs.ctrl, year, stknm) 
    
    n <- array(Obs.stk.nage(biol, ages.error, stk.nage.error, yr+1), dim = c(na, ny+1,it)) # yr+1 in order to be able to calculate F using 'n' ratios
    stock.n(stck)      <- FLQuant(n[,1:ny,], dim = c(na, ny,1,1,1,it), dimnames = dimnames(stock.n(stck)))
    stock.wt(stck)     <- FLQuant(Obs.stk.wgt(biol, ages.error, stk.wgt.error, yr), dim = c(na, ny,1,1,1,it), dimnames = dimnames(stock.n(stck)))
    stock(stck)        <- quantSums(unitSums(seasonSums(stck@stock.n*stck@stock.wt)))
       
    units(stck@harvest) <- 'f'
   
    stck@harvest[-c(na-1,na),] <- log(n[-c(na-1,na),-year,]/n[-c(1,na),-1,]) - stck@m[-c(na-1,na),1:ny,,,,1:it,drop=T]

 #   n. <- array(stck@stock.n[drop=T], dim = c(na,year-1,it))      # [na,ny,it]
 #   m. <- array(stck@m[drop=T], dim = c(na,year-1,it))            # [na,ny,it]
#    c. <- array(stck@catch.n[drop=T], dim = c(na,year-1,it))      # [na,ny,it]
    
     
#    fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
        
#    for(y in 1:ny){
 #        for(a in (na-1):na){
#            for(i in 1:it){
      #      print(i)
#                 zz <- try(ifelse(n.[a,y,i] == 0 | c.[a,y,i] == 0, 0,
 #                                               uniroot(fobj, lower = 0, upper = 1e6, n = n.[a,y,i], m=m.[a,y,i], c = c.[a,y,i])$root))  
#                        stck@harvest[a,y,,,,i] <- ifelse(is.numeric(zz), zz, c(stck@harvest[na-2,y,,,,i]))            }
#        }
#    }

    stck@harvest[c(na-1,na),] <- 0 
    stck@harvest[stck@harvest<0] <- 0

    return(stck)
 }


#-------------------------------------------------------------------------------    
# bio2bioDat(biol, obs.ctrl, year, stknm)  
# Global Biomass Observation of global biomass pop
#-------------------------------------------------------------------------------    
bio2bioDat <- function(biol, fleets, advice, obs.ctrl, year, stknm,...){

    yr                <- year
    stknm             <- stknm
    ny                <- yr-1
    it                <- dim(biol@n)[6]
    
    if (is.null(obs.ctrl$TAC.ovrsht))
      stop(paste("'TAC.ovrsht' array not defined for stock '",stknm,"'",sep=""))
    # TAC.ovrsht  can be numeric with dimension [1]  or FLQuant with dimension [1,dim(biol@n)[2],1,1,1,it]     
    # If TAC.ovrsht is numeric => convert it into an FLQuant. 
    if(is.null(dim(obs.ctrl$TAC.ovrsht))) obs.ctrl$TAC.ovrsht <- FLQuant(obs.ctrl$TAC.ovrsht, dim = c(1,dim(biol@n)[2],1,1,1,it))

    stk.bio.error     <- obs.ctrl$stk.bio.error[,1:ny,drop=F]  
    land.bio.error    <- obs.ctrl$land.bio.error[,1:ny,drop=F] 
    disc.bio.error    <- obs.ctrl$disc.bio.error[,1:ny,drop=F]
    TAC.ovrsht        <- obs.ctrl$TAC.ovrsht[,1:ny,drop=F]
    
    for (e in c('stk.bio.error', 'land.bio.error', 'disc.bio.error')) {
      err <- get(e)
      if (is.null(err))
        stop(paste("'",e,"' array not defined for stock '",stknm,"'",sep=""))
      if(dim(err)[1] != 1)
        stop(paste("'",e,"' array, for stock '",stknm,"', must have dim[1]=1",sep=""))      
      if (sum(err<=0)>0 | sum(is.na(err))>0)
        stop(paste("check values in '",e,"' array for stock '",stknm,"' (required values > 0)",sep=""))
    }
   
    stck              <- propagate(as(biol, "FLStock")[,1:ny], it, fill.iter = TRUE)

    landings(stck)     <- Obs.land.bio(fleets, land.bio.error, yr, stknm)
    landings(stck)     <- FLQuant(ifelse(stck@landings > TAC.ovrsht[1,]*advice$TAC[stknm,1:ny], TAC.ovrsht[1,]*advice$TAC[stknm,1:ny], stck@landings),dim=c(1,ny,1,1,1,it),dimnames=list(age='all', year=dimnames(stck@m)[[2]], unit='unique', season='all', area='unique', iter=1:it))
    discards(stck)     <- Obs.disc.bio(fleets, disc.bio.error, yr, stknm)
    catch(stck)        <- stck@landings + stck@discards

    stck@harvest.spwn[] <- 0 
    stck@m.spwn[]       <- 0 
    stck@range[5]       <- ny
    return(stck)
}

#-------------------------------------------------------------------------------    
# bio2bioPop(biol, obs.ctrl, year, stknm)  
# Global Biomass Observation of global biomass pop
#-------------------------------------------------------------------------------    
bio2bioPop <- function(biol, fleets, advice, obs.ctrl, year, stknm,...){

    stk.bio.error   <- obs.ctrl$stk.bio.error  
    land.bio.error  <- obs.ctrl$land.bio.error 
    disc.bio.error  <- obs.ctrl$disc.bio.error 
    yr              <- year
    stknm           <- stknm
    
    if (is.null(stk.bio.error))
      stop(paste("'stk.bio.error' array not defined for stock '",stknm,"'",sep=""))
    if(dim(stk.bio.error)[1] != 1)
      stop(paste("'stk.bio.error' array, for stock '",stknm,"', must have dim[1]=1",sep=""))      
    if (sum(stk.bio.error<=0)>0 | sum(is.na(stk.bio.error))>0)
      stop(paste("check values in 'stk.bio.error' array for stock '",stknm,"' (required values > 0)",sep=""))
    
    
    stck            <- bio2bioDat(biol, obs.ctrl, yr, stknm)  
    stock(stck)[]      <- Obs.stk.bio(biol, stk.bio.error, yr)
    harvest(stck)[]    <- stck@catch/stck@stock

    return(stck)
}

                        
#-------------------------------------------------------------------------------    
# age2bioDat(biol, fleets, obs.ctrl, year, stknm)
# Age structured to biomass dynamic population.
#-------------------------------------------------------------------------------   
age2bioDat <- function(biol, fleets, advice, obs.ctrl, year, stknm,...){
                         
    yr         <- year
    
    ny <- yr - 1
    it <- dim(biol@n)[6]
    
    if (is.null(obs.ctrl$TAC.ovrsht))
      stop(paste("'TAC.ovrsht' array not defined for stock '",stknm,"'",sep=""))
    # TAC.ovrsht  can be numeric with dimension [1]  or FLQuant with dimension [1,dim(biol@n)[2],1,1,1,it]     
    # If TAC.ovrsht is numeric => convert it into an FLQuant. 
    if(is.null(dim(obs.ctrl$TAC.ovrsht))) obs.ctrl$TAC.ovrsht <- FLQuant(obs.ctrl$TAC.ovrsht, dim = c(1,dim(biol@n)[2],1,1,1,it))

    land.bio.error      <- obs.ctrl$land.bio.error[,1:ny,drop=F] 
    disc.bio.error      <- obs.ctrl$disc.bio.error[,1:ny,drop=F] 
    TAC.ovrsht          <- obs.ctrl$TAC.ovrsht[,1:ny,drop=F]
    
    for (e in c('land.bio.error', 'disc.bio.error')) {
      err <- get(e)
      if (is.null(err))
        stop(paste("'",e,"' array not defined for stock '",stknm,"'",sep=""))
      if(dim(err)[1] != 1)
        stop(paste("'",e,"' array, for stock '",stknm,"', must have dim[1]=1",sep=""))
      if (sum(err<=0)>0 | sum(is.na(err))>0)
        stop(paste("check values in '",e,"' array for stock '",stknm,"' (required values > 0)",sep=""))
      }
             
    biolbio <- setPlusGroupFLBiol(biol,biol@range[1])
    stck                <- propagate(as(biolbio, "FLStock")[,1:ny], it, fill.iter = TRUE) 
    dimnames(stck)      <- list(age="all")

    landings(stck)      <-  FLQuant(Obs.land.bio(fleets, land.bio.error, yr, stknm),dim= c(1,ny,1,1,1,it), dimnames = dimnames(stck@m))
    landings(stck)      <- ifelse(unclass(stck@landings) > TAC.ovrsht[1,]*advice$TAC[stknm,1:ny], TAC.ovrsht[1,]*advice$TAC[stknm,1:ny], stck@landings)
    discards(stck)      <- FLQuant(Obs.disc.bio(fleets, disc.bio.error, yr, stknm),dim= c(1,ny,1,1,1,it), dimnames = dimnames(stck@m))
    catch(stck)         <- stck@landings + stck@discards
    
    cn <- stck@catch 
    ln <- stck@landings
    dn <- stck@discards
    
    names(dimnames(cn))[1] <- names(dimnames(ln))[1] <- names(dimnames(dn))[1] <- 'age'

    catch.n(stck)       <- cn
    landings.n(stck)    <- ln 
    discards.n(stck)    <- dn 
    
    stck@stock <- cn
    stck@harvest <- cn
    stck@stock[] <- NA
    stck@harvest[] <- NA
    
    catch.wt(stck)      <- discards.wt(stck) <- landings.wt(stck) <- 1 
    stck@harvest.spwn[] <- 0 #FLQuant(NA,dim=c(1,ny,1,1,1,it),dimnames=list(age=1, year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))
    stck@m.spwn[]       <- 0 #FLQuant(NA,dim=c(1,ny,1,1,1,it),dimnames=list(age=1, year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))
    stck@mat[]          <- 1
    stck@range[5]       <- stck@range[4]+ny-1
    stck@range[6:7]     <- 0
    return(stck)
}


#-------------------------------------------------------------------------------    
# age2bioPop(biol, fleets, obs.ctrl, year, stknm)
# Age structured to biomass dynamic population.
#-------------------------------------------------------------------------------   
age2bioPop <- function(biol, fleets, advice, obs.ctrl, year, stknm,...){
                         
    stk.bio.error     <- obs.ctrl$stk.bio.error  

    yr                <- year
    
    if (is.null(stk.bio.error))
      stop(paste("'stk.bio.error' array not defined for stock '",stknm,"'",sep=""))
    if(dim(stk.bio.error)[1] != 1)
      stop(paste("'stk.bio.error' array, for stock '",stknm,"', must have dim[1]=1",sep=""))
    if (sum(stk.bio.error<=0)>0 | sum(is.na(stk.bio.error))>0)
      stop(paste("check values in 'stk.bio.error' array for stock '",stknm,"' (required values > 0)",sep=""))
    
         
    stck              <- age2bioDat(biol, fleets, advice, obs.ctrl, yr, stknm)    
    stock(stck)[]       <- Obs.stk.bio(biol, stk.bio.error, yr)
    harvest(stck)[]     <- stck@catch/stck@stock
    
    stck@range[5]     <- yr-1
    return(stck)
}


#-------------------------------------------------------------------------------    
# ageInd(biol, index, obs.ctrl, year, stknm)
# Age structured index
#-------------------------------------------------------------------------------   
# The index is already updated up to year [y-2], we have to update year [y-1].
ageInd <- function(biol, index, obs.ctrl, year, stknm,...){
    
    it <- dim(biol@n)[6]
    ns <- dim(biol@n)[4]
    na <- dim(biol@n)[1]
    ny <- dim(biol@n)[2]
    
    ages.error <- obs.ctrl[['ages.error']]
    
    if(is.null(ages.error)){
        ages.error <- array(0,dim = c(na, na, ny,it), dimnames = list(dimnames(biol@n)[[1]], dimnames(biol@n)[[1]], dimnames(biol@n)[[2]], dimnames(biol@n)[[6]]))
        for(a in 1:na) ages.error[a,a,,] <- 1
    }
    
    if(dim(ages.error)[1] != na | dim(ages.error)[2] != na)
         stop("ages.error array must have dim[1:2] identical to number of ages in stock")
     if(any(round(apply(ages.error,c(1,3:4), sum),2) != 1))
         stop("Some rows in ages.error array  don't add up to 1")
     
    # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
    yrnm   <- dimnames(biol@n)[[2]][year]   
    yrnm.1 <- dimnames(biol@n)[[2]][year-1] 
    
    # season?
    # if the model is seasonal the abundance is taken from the start of the season
    # that corresponds with startf. 
    # sInd: The season from which we are goind to calculate the index. by default sInd = 1.
    sInd <- 1
    if(ns > 1){
        st       <- index@range['startf']
        seasInt  <- seq(0,1,length = ns+1)
        sInd     <- findInterval(st,seasInt)  
    }
    
    ages.ind <- dimnames(index@index.var)[[1]]
    ages.n   <- dimnames(biol@n)[[1]]
    ages.sel <- which(ages.n %in% ages.ind)
    
    for(i in 1:it){
        N <- unitSums(biol@n[ages.sel,yrnm.1,,sInd,,i])
        index@index[,yrnm.1,,,,i] <- (N%*%ages.error[ages.sel,ages.sel,yrnm.1,i])*
                                       index@index.q[,yrnm.1,,,,i, drop=T]*
                                       index@index.var[,yrnm.1,,,,i, drop=T]
    }
    
    return(index)     
}



#-------------------------------------------------------------------------------    
# bioInd(biol, index, obs.ctrl, year, stknm)
# index aggregated in biomass.
#-------------------------------------------------------------------------------   
bioInd <- function(biol, index, obs.ctrl, year, stknm,...){
    
    it <- dim(biol@n)[6]
    ns <- dim(biol@n)[4]
    
    # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
    yrnm   <- dimnames(biol@n)[[2]][year]   
    yrnm.1 <- dimnames(biol@n)[[2]][year-1] 
    
    # season?
    # if the model is seasonal the abundance is taken from the start of the season
    # that corresponds with startf. 
    # sInd: The season from which we are goind to calculate the index. by default sInd = 1.
    sInd <- 1
    if(ns > 1){
        st       <- index@range['startf']
        seasInt  <- seq(0,1,length = ns+1)
        sInd     <- findInterval(st,seasInt)    
    }

    B <- apply(biol@n[,yrnm.1,,sInd,,]*biol@wt[,yrnm.1,,sInd,,],c(2,6),sum)
    index@index[,yrnm.1] <- B*index@index.q[,yrnm.1]*index@index.var[,yrnm.1]
    
    return(index)     
}


#-------------------------------------------------------------------------------    
# SSB index
# ssbInd(biol, index, obs.ctrl, year, stknm)
# index aggregated in biomass.
#-------------------------------------------------------------------------------   
ssbInd <- function(biol, fleets, index, obs.ctrl, year, season,...){
  
  it <- dim(biol@n)[6]
  ns <- dim(biol@n)[4]
  
  # Year  => Character, because the year dimension in indices does not necessarily coincide with year dimension in biol.
  yrnm.1 <- dimnames(biol@n)[[2]][year-1]
  
  # season?
  # Spawning season: determined by stk@m.spwn and stk@harvest.spwn
  sInd <- obs.ctrl$sInd
  # total catch in year [yr] season [ns].
  n.s2 <- biol@n[,yrnm.1,,sInd,]*exp(-biol@m[,yrnm.1,,sInd,])-catchStock(fleets,name(biol))[,yrnm.1,,sInd,]*exp(-biol@m[,yrnm.1,,sInd,]/2)
  fval <- log(biol@n[,yrnm.1,,sInd,]/n.s2) - biol@m[,yrnm.1,,sInd,]
  ssb.stk <- quantSums( biol@n[,yrnm.1,,sInd,]*exp(-(biol@m[,yrnm.1,,sInd,]+fval)*biol@spwn[,yrnm.1,,sInd,])*
                          biol@wt[,yrnm.1,,sInd,]*biol@fec[,yrnm.1,,sInd,])
  
  index@index[,yrnm.1] <- ssb.stk*index@index.q[,yrnm.1]*index@index.var[,yrnm.1]
  
  return(index)     
}


#-------------------------------------------------------------------------------    
# B1,2+ index (in mass)
# cbbmInd(biol, index, obs.ctrl, year, season, stknm)
# index aggregated in biomass.
#-------------------------------------------------------------------------------   
cbbmInd <- function(biol, index, obs.ctrl, year, season,...){
  
  it <- dim(biol@n)[6]
  ns <- dim(biol@n)[4]
  na <- dim(biol@n)[1]
  age1.pos <- 2
  
  # Year  => Character, because the year dimension in indices does not necessarily coincide with year dimension in biol.
  yrnm.1 <- dimnames(biol@n)[[2]][year-1] 
  
  if (season==ns) {
    # season?
    # sInd: The season from which we are goind to calculate the index. 
    # By default sInd = ns. To give B1 at the begging of the following year.
    # IF season = 1 THEN age groups move to the next. 
    sInd <- ns
    # total catch in year [yr] season [ns].
    catch.n <- catchStock(fleets,name(biol))[,yrnm.1,,sInd,]
    # middle ages      # for unit == ss  and age = 1, it will be equal NA but be updated after with SRsim.
    biol.n <- (biol@n[,yrnm.1,,sInd,]*exp(-biol@m[,yrnm.1,,sInd,]/2) - catch.n)*exp(-biol@m[,yrnm.1,,sInd,]/2)
    
    B <- index@index[,yrnm.1,,,]*NA
    B[1,] <- biol.n[age1.pos-1,]*obs.ctrl$wage['age1',]
    B[2,] <- quantSums(biol.n[(age1.pos):na,])*obs.ctrl$wage['age2plus',]
    index@index[,yrnm.1,] <- B*index@index.q[,yrnm.1,]*index@index.var[,yrnm.1,]
  }
  
  return(index)     
}
