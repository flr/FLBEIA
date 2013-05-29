# author: Ruben Roa (2010/11/XY to 2011/03/29)
# changed: Dorleta Garcia (2011/04/26)

################################################################################
##################### AGE-STRUCTURED POP OBSERVED BY AGE #######################
################################################################################
## OPERATING ON OBJECTS OF CLASS FLBiols
#Age reading observation error (observation) and total N observation error (estimation)
#error.ages : an array of probabilities of age assignment of dim =c(na, na, ny, iter)
#varia.ntot : an array of multiplicative errors (total abundance estimation error)
#yr : integer, the year of assessment/management
#Details
#error.ages is an input array that for each yr and iter provides a square matrix whose column elements quantify the probability 
#that a fish of age 'a' be classified as having any age in the seq min(age):max(age). For example a column such as c(0,0.25,0.5,0.25,0,0) 
#for the third age class sets the probability of correct age classification at 0.5 whereas the probability that a fish of age 'a' be 
#classified in the 1st age class is zero, second age class 0.25, and so on. A change of the na*na matrix over the years would mean
#improvement or worsening of the age reading error. A change of the matrix over the iterations would represent uncertainty about the errors
#of age assignment, error in knowing the error if you will.
#The na*na matrix in error.ages is conditioned on column totals (which should add up to 1) thus it is age-assignment error exclusively.
#If the na*na matrix is the identity matrix, age assignment error is suppressed.
#If the array varia.ntot has errors such that they induce negative numbers in the observation, 
#then the function will stop with an error message.
#If varia.ntot = 0 total abundance estimation error is suppressed.
Obs.ages <- function(biol, error.ages, varia.ntot, yr)
            {                                                    
            na  <- dim(biol@n)[1]
            ny  <- yr  -1
            ns <- dim(biol@n)[4]
            it  <- dim(biol@n)[6]
            bio.num <- unitSums(biol@n[,,,1,,]) #sum thru units of numbers in first season
            
            bio.num[1,,,,,] <-  0
            for(ss in 1:ns) bio.num[1,,,,,] <- bio.num[1,,,,,] + biol@n[1,,ss,ss,,] # unitSums(seasonSums(biol@n[1,,,,,])) #replace recruitment by sum of recruitment thru time steps
           
            bio.num <- bio.num[,1:ny,,,,]
 
            obs.num <- bio.num; obs.num[] <- 0
            for(j in 1:ny){
              for(n in 1:it){
                obs.num[,j,,,,n] <-bio.num[,j,,,,n]%*%error.ages[,,j,n]
                }
            }

            obs.num            <- obs.num*varia.ntot[,1:ny]
             return(obs.num)
            }
#Natural mortality random variation
#varia.mort : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
#yr : integer, the year the stock is observed from
#Details: If varia.mort is 0 then there is no uncertainty in M.
Var.mort <- function(biol,  error.ages, varia.mort,  yr)
             {
    ny   <- yr  -1
    mort.real <- seasonMeans(unitMeans(biol@m))[,1:ny,]
    mort.obs  <- mort.real
    it <- dim(mort.real)[6]
             
    for(i in 1:it){
        for(y in 1:ny){
            mort.obs[,y,,,,i] <- mort.real[,y,,,,i]%*%error.ages[,,y,i] # Here we should weight byu the real age composition of the population??? 
        }
    }
    
    mort.obs  <- mort.obs*varia.mort[,1:ny]
    return(mort.obs)
}

#Weight at age observation error      IN TTHE POPULATION
#var.mwgt : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
#yr : integer, the year the stock is observed from
#Details
#If var.mwgt is 1 then there is no observation error in mean weight at age
Obs.mwgt <- function(biol, error.ages, varia.mwgt, yr){
    ny   <- yr  -1
    mwgt.real  <- seasonMeans(unitMeans(biol@wt))[,1:ny]
    mwgt.obs  <- mwgt.real
    it <- dim(mwgt.real)[6]
             
    for(i in 1:it){                                               # Here we should weight byu the real age composition of the population??? 
        for(y in 1:ny){
            mwgt.obs[,y,,,,i] <- mwgt.real[,y,,,,i]%*%error.ages[,,y,i]
        }
    }
    mwgt.obs <- mwgt.obs*varia.mwgt[,1:ny]
    return(mwgt.obs)
}

#Fec at age observation error
#yr : integer, the year the stock is observed from
#Details
#Observation error of proportions are not multiplicative nor multiplicative. Thus to introduce observation error in proportion
#mature it is necessary to replace the real values in totto with the observed values.
#In creating the replacement observed values the mean will be given by the real values and the shape2 parameter of the Beta distribution.
#bot.age and top.age allows to constrain the uncertainty to a range of intermediate ages, taking advantage of the fact that very young fish 
#(i.e. age 0) will not be reproducing no matter what, and obviously adult fish will be sexually mature (ceteris paribus).
#If bot.age and top.age are equal then maturity will go from 0 at age (bot.age - 1) to 1 at age (top.age)
#varia.mwgt : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)

Obs.fec  <- function(biol, error.ages, varia.fec, yr){
    
    ny   <- yr -1
    
    fec.real  <- seasonMeans(unitMeans(biol@fec))[,1:ny]
    fec.obs  <- fec.real
    it <- dim(fec.real)[6]
    
    for(i in 1:it){
        for(y in 1:ny){
            fec.obs[,y,,,,i] <- fec.real[,y,,,,i]%*%error.ages[,,y,i] 
        }
    }
    fec.obs  <- fec.obs*varia.fec[,1:ny]
    
    return(fec.obs)
}
#Total numbers per year observation error  - Global models
#varia.nyr : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
#yr : integer, the year the stock is observed from
#Details: If varia.nyr is 0 then there is no uncertainty in total numbers per year.
Var.nyr <- function(biol, varia.nyr, yr){
    ny   <- yr -1
    nyr  <- biol@n[,,,1,,]
    nyr  <- nyr*varia.nyr[,1:ny]
    return(nyr)
}
##OPERATING ON OBJECTS OF CLASS FLFleetExt
#Landings in numbers at age observation error
#fleets : an object of class FLFleetsExt
#error.ages : an array of probabilities of age assignment of dim =c(na, na, ny, iter)
#varia.ltot : an array of multiplicative errors 
#lrep.bias : a positive number (equal for all years and ages) or a vector of such numbers f (change over the years) or an array
#(change over the years and ages) of multiplicative systematic under- or over-reporting.
#yr : integer, the year the stock is observed from
#stknm : character, the name of the stock
#Details
#error.ages is an input array that for each yr and iter provides a square matrix whose column elements quantify the probability 
#that a fish of age 'a' be classified as having any age in the seq min(age):max(age). For example a column such as c(0,0.25,0.5,0.25,0,0) 
#for the third age class sets the probability of correct age classification at 0.5 whereas the probability that a fish of age 'a' be 
#classified in the 1st age class is zero, second age class 0.25, and so on. A change of the na*na matrix over the years would mean
#improvement or worsening of the age reading error. A change of the matrix over the iterations would represent uncertainty about the errors
#of age assignment, error in knowing the error if you will.
#The na*na matrix in error.ages is conditioned on column totals (which should add up to 1) thus it is age-assignment error exclusively.
#If the na*na matrix is the identity matrix, age assignment error is suppressed.
#If the array varia.ntot has errors such that they induce negative numbers in the observation, 
#then the function will stop with an error message.
#If varia.ntot = 0 total abundance estimation error is suppressed.
Obs.laage <- function(fleets, error.ages, varia.ltot, obs.lwta, yr, stknm){

    ny     <- yr -1

    la.num <- unitSums(seasonSums(landStock(fleets,stknm)))[,1:ny] #sum thru units of numbers of all season

    na     <- dim(la.num)[1]
    it     <- dim(la.num)[6]     
       
    ola.num <- array(0,c(na,ny,1,1,1,it))
    mwgt.real <- seasonMeans(unitMeans(wtalStock(fleets,stknm)))[,1:ny]

    for(j in 1:ny){
       for(n in 1:it){
           ola.num[,j,,,,n]  <- la.num[,j,,,,n]%*%error.ages[,,j,n]
        }
    }
     
    # Correct the landings at age so that _new_ total landings are equal to 
    # _real_ total landings (corrected by dga)
    Lnew  <- apply(ola.num*obs.lwta, c(2:6), sum)  # [ny,1,1,1,it] (the factors of the * are arrays)
    Lreal <- apply(la.num*mwgt.real, c(2:6), sum)  # [1,ny,1,1,1,it]  (the factors of the * are FLQ-s)
    Lrat <- Lreal/array(Lnew, dim = c(1, dim(Lnew)))
    
    Lrat[is.na(Lrat)] <- 0
    
    ola.num <- sweep(ola.num, 2:6, Lrat, "*")
    
    # Apply the error in the observation of the landings due to misreporting.
    ola.num          <- ola.num*varia.ltot[,1:ny]
    
    return(ola.num)
}
#Weight at age in landings observation error
#fleets : an object of class FLFleetExt
#varia.mwgt : an array of dim = c(na,ny,1,1,1,it)
#yr : integer, last year of observation
#stknm : character, name of the stock
Obs.wtal <- function(fleets, error.ages, varia.lwgt, laa.obs, yr, stknm){

     ny   <- yr -1
     laa <- unitSums(seasonSums(landStock(fleets,stknm)))[,1:ny]
     
     it <- dim(fleets[[1]]@effort)[6]
     
     mwgt.real <- seasonMeans(unitMeans(wtalStock(fleets,stknm)))[,1:ny]
     mwgt.obs  <- mwgt.real
     
     for(i in 1:it){
        for(y in 1:ny){
       # [rr]     mwgt.obs[,y,,,,i] <- ((mwgt.real[,y,,,,i,drop=T]*laa[,y,,,,i,drop=T])%*%error.ages[,,y,i])/laa.obs[,y,,,,i, drop=T]
            mwgt.obs[,y,,,,i] <- mwgt.real[,y,,,,i,drop=T]%*%error.ages[,,y,i]

        }
     }
    
     mwgt.obs <- mwgt.obs*varia.lwgt[,1:ny]
     return(mwgt.obs)
}
#Discards in numbers at age observation error
#fleets : an object of class FLFleetsExt
#error.ages : an array of probabilities of age assignment of dim =c(na, na, ny, iter)
#varia.ntot : an array of multiplicative errors (total abundance estimation error)
#drep.bias : a positive number or a vector of such numbers for multiplicative systematic under- or over-reporting
#yr : integer, the year the stock is observed from
#stknm : character, the name of the stock
#Details
#error.ages is an input array that for each yr and iter provides a square matrix whose column elements quantify the probability 
#that a fish of age 'a' be classified as having any age in the seq min(age):max(age). For example a column such as c(0,0.25,0.5,0.25,0,0) 
#for the third age class sets the probability of correct age classification at 0.5 whereas the probability that a fish of age 'a' be 
#classified in the 1st age class is zero, second age class 0.25, and so on. A change of the na*na matrix over the years would mean
#improvement or worsening of the age reading error. A change of the matrix over the iterations would represent uncertainty about the errors
#of age assignment, error in knowing the error if you will.
#The na*na matrix in error.ages is conditioned on column totals (which should add up to 1) thus it is age-assignment error exclusively.
#If the na*na matrix is the identity matrix, age assignment error is suppressed.
#If the array varia.ntot has errors such that they induce negative numbers in the observation, 
#then the function will stop with an error message.
#If varia.ntot = 0 total abundance estimation error is suppressed.
Obs.daage <- function(fleets, error.ages, varia.dtot, obs.dwta, yr, stknm){
    #sum thru units of numbers of all season            
    
    ny     <- yr -1 
    da.num <- unitSums(seasonSums(discnStock(fleets,stknm)))[,1:ny]
    it     <- dim(da.num)[6]
    na     <- dim(da.num)[1]
    
    oda.num <- da.num; oda.num[] <- 0
    
    mwgt.real <- seasonMeans(unitMeans(wtadStock(fleets,stknm)))[,1:ny]

    
    for(j in 1:ny){
        for(n in 1:it){
            oda.num[,j,,,,n] <- da.num[,j,,,,n]%*%error.ages[,,j,n]
        }
    }
    
         
    # Correct the discards at age so that _new_ total discards are equal to 
    # _real_ total discards (corrected by dga)
    Dnew  <- apply(oda.num*obs.dwta, c(2:6), sum)  # [1,ny,1,1,1,it]
    Dreal <- apply(da.num*mwgt.real, c(2:6), sum)  # [1,ny,1,1,1,it] 
    Drat <- Dreal/Dnew

    Drat[is.na(Drat),] <- 1  # Drat == NA => is because Dnew = Dreal == 0.
    
    oda.num <- sweep(oda.num, 2:6, Drat, "*")
    
    # Apply the error in the observation of the discards due to misreporting.
    oda.num          <- oda.num*varia.dtot[,1:ny]
    
    return(oda.num)
}
#Weight at age in discards observation error
#fleets : an object of class FLFleetExt
#varia.mwgt : an array of dim = c(na,ny,1,1,1,it)
#yr: integer, last year of observation
#stknm : character, name of the stock
Obs.wtad <- function(fleets, error.ages, varia.dwgt,  daa.obs, yr, stknm){
    
    ny   <- yr -1
    daa <- unitSums(seasonSums(discStock(fleets,stknm)))[,1:ny]
    
    it <- dim(fleets[[1]]@effort)[6]
    
    mwgt.real <- seasonMeans(unitMeans(wtadStock(fleets,stknm)))[,1:ny]
    mwgt.obs  <- mwgt.real
    
    for(i in 1:it){
       for(y in 1:ny){
      # [rr]     mwgt.obs[,y,,,,i] <- ((mwgt.real[,y,,,,i]*daa[,y,,,,i])%*%error.ages[,,y,i])/daa.obs[,y,,,,i,drop=T]
            mwgt.obs[,y,,,,i] <- mwgt.real[,y,,,,i]%*%error.ages[,,y,i]
       }
    }
    
    mwgt.obs <- mwgt.obs*varia.dwgt[,1:ny]
    
    return(mwgt.obs)
}
#Harvest at age from solving Baranov eq. with uniroot
Obs.FSolver <- function(Fay,May,Nay,Cay)
       {
       Cay-(Fay/(Fay+May))*(1-exp(-Fay-May))*Nay
       }

################################################################################
################### GLOBAL POP OBSERVED BY GLOBAL BIOMASS ######################
################################################################################
## OPERATING ON OBJECTS OF CLASS FLBiols
#Total biomass observation error
#biol : an object of class FLBiols
#varia.btot <- a number or a numeric vector of coefficients for systematic and/or random error in the observation of total biomass
#yr : integer, the year the stock is observed from
Obs.btot <- function(biol, varia.btot, yr){
     btot     <- unitSums(biol@n[,,,1,,]*biol@wt[,,,1,,])
     ny       <- yr - 1 
     btot     <- btot*varia.btot[,1:ny]
     return(btot)
}
##OPERATING ON OBJECTS OF CLASS FLFleetExt
#Total landings observation error
#fleets : an object of class FLFleetsExt
#varia.tland : a number or a numeric vector of multiplicative deviates for random error in the observation of landings
#lrep.bias : a positive number or a vector of such numbers for multiplicative systematic under- or over-reporting
#yr : integer, the year the stock is observed from
#stknm : character, identifier of the stock
#Details
#A value of 0 for varia.tland and a value of 1 for lrep.bias means that the observation of the total landings is perfect.
Obs.tland <- function(fleets, varia.tland,  yr, stknm){                                                                       
    ny        <- yr - 1
    tland     <- unitSums(seasonSums(tlandStock(fleets, stknm)))[,1:ny]
    tland     <- tland*varia.tland[,1:ny]

    return(tland)
}
#Total discards observation error
#fleets : an object of class FLFleetsExt
#varia.tdisc : a number or a numeric vector of multiplicative deviates for random error in the observation of landings
#drep.bias : a positive number or a vector of such numbers for multiplicative systematic under- or over-reporting
#yr : integer, the year the stock is observed from
#stknm : character, identifier of the stock
#Details
#A value of 0 for varia.tdisc and a value of 1 for drep.bias means that the observation of the total discards is perfect.
Obs.tdisc <- function(fleets, varia.tdisc, yr, stknm){
     ny        <- yr -1
     tdisc     <- unitSums(seasonSums(tdiscStock(fleets, stknm)))[,1:ny]
     tdisc     <- tdisc*varia.tdisc[,1:ny]

     return(tdisc)
}

################################################################################
############# AGE-STRUCTURED POP OBSERVED BY GLOBAL BIOMASS ####################
################################################################################

#OPERATING ON OBJECTS OF CLASS FLFleetExt
#Total landings observation error
#fleets : an object of class FLFleetsExt
#varia.tland : a number or a numeric vector of multiplicative deviates for random error in the observation of landings
#lrep.bias : a positive number or a vector of such numbers for multiplicative systematic under- or over-reporting
#yr : integer, the year the stock is observed from
#stknm : character, identifier of the stock
#Details
#A value of 0 for varia.tland and a value of 1 for lrep.bias means that the observation of the total landings is perfect.
Obs.tlaas <- function(fleets, varia.tland,  yr, stknm){
    ny        <- yr-1 
    tlaas     <- unitSums(seasonSums(tlandStock(fleets, stknm)))[,1:ny]
    tlaas     <- tlaas*varia.tland[,1:ny]
    return(tlaas)
}
#Total discards observation error
#fleets : an object of class FLFleetsExt
#varia.tdias : a number or a numeric vector of multiplicative deviates for random error in the observation of landings
#drep.bias : a positive number or a vector of such numbers for multiplicative systematic under- or over-reporting
#yr : integer, the year the stock is observed from
#stknm : character, identifier of the stock
#Details
#A value of 0 for varia.tdisc and a value of 1 for drep.bias means that the observation of the total discards is perfect.
Obs.tdias <- function(fleets, varia.tdisc,  yr, stknm){
    ny        <- yr - 1
    tdias     <- unitSums(seasonSums(tdiscStock(fleets, stknm)))[,1:ny]
    tdias     <- tdias*varia.tdisc[,1:ny]

    return(tdias)
    }
#END