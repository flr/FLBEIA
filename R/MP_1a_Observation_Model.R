#-------------------------------------------------------------------------------
#                       OBSERVATION MODEL FUNCTIONS
#
#   - perfectObservation: Variables in FLStock are observed without error.  
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
#
# Dorleta GarcYYYa
# Created: 09/12/2010 14:37:43
# Changed: 15/04/2013 10:59:56 (correct a bug in PerfectObs problem with seasons) 
#(rroa's functions inserted)
#         2013-06-07 12:20:04 Sonia Sanchez - code revised and names changed for coherence
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# perfectObservation(biol, fleets, covars, obs.ctrl, year = 1, season = 1)
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
    res <- as(biol, 'FLStock')[,1:(year-1),1,1]
        
    res@range[c(1:3,6:7)] <- biol@range[c(1:3,6:7)]
    names(res@range[6:7]) <- c('minfbar', 'maxfbar')
        
    res@discards.wt <- res@landings.wt <- res@catch.wt <- res@stock.wt
        
    # "stock.n":  FIRST SEASON and SUM ALONG UNITS except recruitment
    # rec = n[1,,1,1] + n[1,,2,2] + n[1,,3,3] + n[1,,4,4]
    # n up to (year) to use it after in the 'f' calculation.
    n <- unitSums(biol@n)[,1:year,,1]
        
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
         
    res@stock.n[] <- n[,1:(year-1)]
    
    res@stock[] <- quantSums(res@stock.n*res@stock.wt)
        
    # SUM ALONG SEASONS AND FIRST UNIT: "m"
    res@m[]      <- seasonSums(biol@m)[,1:(year-1),1,]
    res@m.spwn[] <- seasonSums(biol@spwn)[,1:(year-1),1,]/ns
        
    # SUM ALONG UNITS AND SEASONS, OBTAINED FROM FLFLEETS: 
    # "catch", "catch.n", "discards"     "discards.n" "landings"     "landings.n"
    res@landings.n[] <- apply(landStock(fleets, st), c(1:2,6),sum)[,1:(year-1),]
    res@discards.n[] <- apply(discStock(fleets, st), c(1:2,6),sum)[,1:(year-1),]
    res@catch.n[]    <- res@discards.n + res@landings.n
    res@landings[]   <- quantSums(res@landings.n*res@landings.wt)
    res@discards[]   <- quantSums(res@discards.n*res@discards.wt)
    res@catch[]      <- quantSums(res@catch.n*res@catch.wt)
        
    # harvest: * if age structured calculate it from 'n'.
    #          * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
    if(na == 1){
        res@harvest[] <- res@catch/(res@stock.n*res@stock.wt)
        units(res@harvest) <- 'hr'
    } else{
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
                    n.[a,y,i] - c.[a,y,i]
                    if(n.[a,y,i] < c.[a,y,i]) res@harvest[a,y,,,,i] <- 10
                    
                    else{
                        zz <- try(ifelse(n.[a,y,i] == 0 | c.[a,y,i] == 0, 0,
                                                uniroot(fobj, lower = 0, upper = 10, n = n.[a,y,i], m=m.[a,y,i], c = c.[a,y,i])$root))  
                        res@harvest[a,y,,,,i] <- ifelse(is.numeric(zz), zz, res@harvest[ai-1,y,,,,i] )
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
        for(s in c(1:ns)) c.perc[,,,s,] <- ctot.age[,,,s,]/ctot
        biol.spwn <- unitMeans(biol@spwn[,1:(year-1),,])
        res@harvest.spwn[] <- seasonSums(c.perc[,1:(year-1),,]*biol.spwn)

    }

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
                                                
    # TAC.ovrsht  can be numeric with dimension [1]  or FLQuant with dimension [1,dim(biol@n)[2],1,1,1,it]     
    # If TAC.ovrsht is numeric => convert it into an FLQuant. 
    if(is.null(dim(obs.ctrl$TAC.ovrsht))) obs.ctrl$TAC.ovrsht <- FLQuant(obs.ctrl$TAC.ovrsht, dim = c(1,dim(biol@n)[2],1,1,1,it))
                                                                
    ages.error        <- obs.ctrl$ages.error
    nmort.error       <- obs.ctrl$nmort.error[,1:ny]
    fec.error         <- obs.ctrl$fec.error[,1:ny]
    land.wgt.error    <- obs.ctrl$land.wgt.error[,1:ny]
    disc.wgt.error    <- obs.ctrl$disc.wgt.error[,1:ny]
    land.nage.error   <- obs.ctrl$land.nage.error[,1:ny]
    disc.nage.error   <- obs.ctrl$disc.nage.error[,1:ny]   
    TAC.ovrsht        <- obs.ctrl$TAC.ovrsht[,1:ny]
        
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
         
    stck              <- as(biol, "FLStock")[,1:ny,1,1]    

    stck@landings.wt  <- Obs.land.wgt(fleets, ages.error, land.wgt.error, yr, stknm)
    stck@landings.n[] <- Obs.land.nage(fleets, ages.error, land.nage.error, stck@landings.wt, yr, stknm)
    stck@landings     <- quantSums(seasonSums(stck@landings.n*stck@landings.wt))
    
    # In landings.wt the error due to age depends on landings.n, but that on m and mat only depends on error itself
    # because it is suppose that the biological sampling is independent.
    
    stck@m            <- Obs.nmort(biol, ages.error, nmort.error, yr)
    stck@mat          <- Obs.fec(biol, ages.error, fec.error, yr)

    # compare the landings with the advice and depending on TAC.ovrsht report landings. the misresporting is 
    # reported homogeneously.
    stck@landings     <- FLQuant(ifelse(unclass(stck@landings) > TAC.ovrsht*advice$TAC[stknm,1:ny], TAC.ovrsht*advice$TAC[stknm,1:ny], stck@landings))
    
    ovrsht.red        <- stck@landings/quantSums(unitSums(seasonSums(stck@landings.n*stck@landings.wt)))  # [1,ny,,,,it]
    
    ovrsht.red[is.na(ovrsht.red)] <- 1
     
    stck@landings.n   <- sweep(stck@landings.n, 2:6, ovrsht.red, "*") #distributing the overshoot subreporting of bulk landings in biomass equally over ages
   
    stck@discards.wt  <- Obs.disc.wgt(fleets, ages.error, disc.wgt.error,  yr, stknm)
    stck@discards.wt[is.na(stck@discards.wt)] <- stck@landings.wt[is.na(stck@discards.wt)]
    stck@discards.n   <- Obs.disc.nage(fleets, ages.error, disc.nage.error, stck@discards.wt, yr, stknm)
    stck@discards     <- quantSums(seasonSums(stck@discards.n*stck@discards.wt))
    
    stck@catch        <- stck@landings + stck@discards
    stck@catch.n      <- stck@landings.n + stck@discards.n
    stck@catch.wt     <- (stck@landings.n*stck@landings.wt + stck@discards.n*stck@discards.wt)/(stck@landings.n + stck@discards.n)
    
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
    
    n <- Obs.stk.nage(biol, ages.error, stk.nage.error, yr)
    stck@stock.n      <- n[,1:ny,]
    stck@stock.wt     <- Obs.stk.wgt(biol, ages.error, stk.wgt.error, yr)
    stck@stock        <- quantSums(unitSums(seasonSums(stck@stock.n*stck@stock.wt)))
       
    units(stck@harvest) <- 'f'
   
    stck@harvest[-c(na-1,na),] <- log(n[-c(na-1,na),-year]/n[-c(1,na),-1]) - stck@m[-c(na-1,na),]

    n. <- array(stck@stock.n[drop=T], dim = c(na,year-1,it))      # [na,ny,it]
    m. <- array(stck@m[drop=T], dim = c(na,year-1,it))            # [na,ny,it]
    c. <- array(stck@catch.n[drop=T], dim = c(na,year-1,it))      # [na,ny,it]
        
    fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
        
    for(y in 1:ny){
        for(a in (na-1):na){
            for(i in 1:it){
      #      print(i)
                 zz <- try(ifelse(n.[a,y,i] == 0 | c.[a,y,i] == 0, 0,
                                                uniroot(fobj, lower = 0, upper = 10, n = n.[a,y,i], m=m.[a,y,i], c = c.[a,y,i])$root))  
                        stck@harvest[a,y,,,,i] <- ifelse(is.numeric(zz), zz, c(stck@harvest[na-2,y,,,,i]))            }
        }
    }
    
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
    
    # TAC.ovrsht  can be numeric with dimension [1]  or FLQuant with dimension [1,dim(biol@n)[2],1,1,1,it]     
    # If TAC.ovrsht is numeric => convert it into an FLQuant. 
    if(is.null(dim(obs.ctrl$TAC.ovrsht))) obs.ctrl$TAC.ovrsht <- FLQuant(obs.ctrl$TAC.ovrsht, dim = c(1,dim(biol@n)[2],1,1,1,it))

    stk.bio.error     <- obs.ctrl$stk.bio.error[,1:ny]  
    land.bio.error    <- obs.ctrl$land.bio.error[,1:ny] 
    disc.bio.error    <- obs.ctrl$disc.bio.error[,1:ny]
    TAC.ovrsht        <- obs.ctrl$TAC.ovrsht[,1:ny]
    
    for (e in c('stk.bio.error', 'land.bio.error', 'disc.bio.error')) {
      err <- get(e)
      if (is.null(err))
        stop(paste("'",e,"' array not defined for stock '",stknm,"'",sep=""))
      if(dim(err)[1] != 1)
        stop(paste("'",e,"' array, for stock '",stknm,"', must have dim[1]=1",sep=""))      
      if (sum(err<=0)>0 | sum(is.na(err))>0)
        stop(paste("check values in '",e,"' array for stock '",stknm,"' (required values > 0)",sep=""))
    }
   
    stck              <- as(biol, "FLStock")[,1:ny]

    stck@landings     <- Obs.land.bio(fleets, land.bio.error, yr, stknm)
    stck@landings     <- FLQuant(ifelse(stck@landings > TAC.ovrsht*advice$TAC[stknm,1:ny], TAC.ovrsht*advice$TAC[stknm,1:ny], stck@landings),dim=c(1,ny,1,1,1,it),dimnames=list(age='all', year=dimnames(stck@m)[[2]], unit='unique', season='all', area='unique', iter=1:it))
    stck@discards     <- Obs.disc.bio(fleets, disc.bio.error, yr, stknm)
    stck@catch        <- stck@landings + stck@discards

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
    stck@stock      <- Obs.stk.bio(biol, stk.bio.error, yr)
    stck@harvest    <- stck@catch/stck@stock

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
    
    # TAC.ovrsht  can be numeric with dimension [1]  or FLQuant with dimension [1,dim(biol@n)[2],1,1,1,it]     
    # If TAC.ovrsht is numeric => convert it into an FLQuant. 
    if(is.null(dim(obs.ctrl$TAC.ovrsht))) obs.ctrl$TAC.ovrsht <- FLQuant(obs.ctrl$TAC.ovrsht, dim = c(1,dim(biol@n)[2],1,1,1,it))

    land.bio.error      <- obs.ctrl$land.bio.error[,1:ny] 
    disc.bio.error      <- obs.ctrl$disc.bio.error[,1:ny] 
    TAC.ovrsht          <- obs.ctrl$TAC.ovrsht[,1:ny]
    
    for (e in c('land.bio.error', 'disc.bio.error')) {
      err <- get(e)
      if (is.null(err))
        stop(paste("'",e,"' array not defined for stock '",stknm,"'",sep=""))
      if(dim(err)[1] != 1)
        stop(paste("'",e,"' array, for stock '",stknm,"', must have dim[1]=1",sep=""))
      if (sum(err<=0)>0 | sum(is.na(err))>0)
        stop(paste("check values in '",e,"' array for stock '",stknm,"' (required values > 0)",sep=""))
      }
             
    biolbio <- setPlusGroup(biol,biol@range[1])
    stck                <- as(biolbio, "FLStock")[,1:ny]
    stck@landings       <- Obs.land.bio(fleets, land.bio.error, yr, stknm)
    stck@landings[]     <- ifelse(unclass(stck@landings) > TAC.ovrsht*advice$TAC[stknm,1:ny], TAC.ovrsht*advice$TAC[stknm,1:ny], stck@landings)
    stck@discards       <- Obs.disc.bio(fleets, disc.bio.error, yr, stknm)
    stck@catch          <- stck@landings + stck@discards
    stck@catch.n[]      <- stck@catch 
    stck@landings.n[]   <- stck@landings 
    stck@discards.n[]   <- stck@discards 
    stck@catch.wt[]     <- stck@discards.wt[] <- stck@landings.wt[] <- 1 
    stck@harvest.spwn[] <- 0 #FLQuant(NA,dim=c(1,ny,1,1,1,it),dimnames=list(age=1, year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))
    stck@m.spwn[]       <- 0 #FLQuant(NA,dim=c(1,ny,1,1,1,it),dimnames=list(age=1, year=biol@range[4]:ny, unit='unique', season='all', area='unique', iter=1:it))
    stck@mat[]          <- 1
    stck@range[5]       <- ny
    stck@range[6:7]     <- 1
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
    stck@stock        <- Obs.btot(biol, stk.bio.error, yr)
    stck@harvest      <- stck@catch/stck@stock
    
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
