# author: Ruben Roa (2010/11/XY to 2011/03/29)
# changed: Dorleta Garcia (2011/04/26)
#         2013-06-07 12:20:04 Sonia Sanchez - code revised and names changed for coherence

#-------------------------------------------------------------------------------
#               POPULATION OBSERVATION FUNCTIONS
#
# Obs.stk.nage  : Age reading observation error (observation) and total N observation error (estimation)
# Obs.nmort     : Natural mortality random variation
# Obs.stk.wgt   : Weight at age in the stock observation error
# Obs.fec       : Fecundity at age observation error
# Obs.land.nage : Landings in numbers at age observation error
# Obs.land.wgt  : Weight at age in landings observation error
# Obs.disc.nage : Discards in numbers at age observation error
# Obs.disc.wgt  : Weight at age in discards observation error
# Obs.FSolver   : Harvest at age from solving Baranov eq. with uniroot
# Obs.stk.bio   : Total biomass observation error
# Obs.land.bio  : Total landings observation error
# Obs.disc.bio  : Total discards observation error
#-------------------------------------------------------------------------------


################################################################################
##################### AGE-STRUCTURED POP OBSERVED BY AGE #######################
################################################################################

## OPERATING ON OBJECTS OF CLASS FLBiols
##########################################

# Obs.stk.nage
#---------------
# Function to simulate the population observation in numbers at age at the beginning of the year, taking into account the age misclassification
# - Age reading observation error (observation) and total N observation error (estimation)
#Input
# biol          : an an object of class FLBiol
# ages.error    : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# stk.nage.error: an array of multiplicative errors (total abundance estimation error), dim = c(na,ny)
# yr            : integer, the year of assessment/management
#Details
# ages.error is an input array that for each yr and iter provides a square matrix whose column elements quantify the probability 
#that a fish of age 'a' is classified as having any age in the seq min(age):max(age). For example a column such as c(0,0.25,0.5,0.25,0,0) 
#for the third age class sets the probability of correct age classification at 0.5 whereas the probability that a fish of age 'a' be 
#classified in the 1st age class is zero, second age class 0.25, and so on. A change of the na*na matrix over the years would mean
#improvement or worsening of the age reading error. A change of the matrix over the iterations would represent uncertainty about the errors
#of age assignment, error in knowing the error if you will.
# The na*na matrix in ages.error is conditioned on column totals (which should add up to 1) thus it is age-assignment error exclusively.
# If the na*na matrix is the identity matrix, age assignment error is suppressed.
# If stk.nage.error = 1 numbers at age estimation error is suppressed.
Obs.stk.nage <- function(biol, ages.error, stk.nage.error, yr)
            {                                                    
            na  <- dim(biol@n)[1]
            ny  <- yr  -1
            ns <- dim(biol@n)[4]
            it  <- dim(biol@n)[6]
            bio.num <- unitSums(biol@n[,,,1,,]) #sum thru units of numbers in first season
            
            bio.num[1,,,,,] <-  0
            for(ss in 1:ns) bio.num[1,,,,,] <- bio.num[1,,,,,] + biol@n[1,,ss,ss,,] #replace recruitment by sum of recruitment thru time steps
           
            bio.num <- bio.num[,1:ny,,,,]
 
            obs.num <- bio.num; obs.num[] <- 0
            for(j in 1:ny){
              for(n in 1:it){
                obs.num[,j,,,,n] <-bio.num[,j,,,,n]%*%ages.error[,,j,n]
                }
            }

            obs.num            <- obs.num*stk.nage.error[,1:ny]
             return(obs.num)
            }

# Obs.nmort
#---------------
# Function to simulate the observation of the mean weight at age in population (seasonal mean), taking into account the age misclassification
# - Natural mortality random variation
#Input
# biol        : an an object of class FLBiol
# ages.error  : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# nmort.error : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# yr          : integer, the year the stock is observed from
#Details
# If nmort.error is 1 then there is no uncertainty in M.
Obs.nmort <- function(biol,  ages.error, nmort.error,  yr)
             {
    ny   <- yr  -1
    mort.real <- seasonSums(unitMeans(biol@m))[,1:ny,]
    mort.obs  <- mort.real
    it <- dim(mort.real)[6]
             
    for(i in 1:it){
        for(y in 1:ny){
            mort.obs[,y,,,,i] <- mort.real[,y,,,,i]%*%ages.error[,,y,i] #! Here we should weight by the real age composition of the population??? 
        }
    }
    
    mort.obs  <- mort.obs*nmort.error[,1:ny]
    return(mort.obs)
}

# Obs.stk.wgt
#---------------
# Function to simulate the observation of the mean weight at age in population (seasonal mean), taking into account the age misclassification
# - Weight at age in the stock observation error
#Input
# biol          : an an object of class FLBiol
# ages.error    : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# stk.wgt.error : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# yr            : integer, the year the stock is observed from
#Details
# If stk.wgt.error is 1 then there is no observation error in mean weight at age in the stock.
Obs.stk.wgt <- function(biol, ages.error, stk.wgt.error, yr){
    ny   <- yr  -1
    mwgt.real  <- seasonMeans(unitMeans(biol@wt))[,1:ny]
    mwgt.obs  <- mwgt.real
    it <- dim(mwgt.real)[6]
             
    for(i in 1:it){                                               #! Here we should weight by the real age composition of the population??? 
        for(y in 1:ny){
            mwgt.obs[,y,,,,i] <- mwgt.real[,y,,,,i]%*%ages.error[,,y,i]
        }
    }
    mwgt.obs <- mwgt.obs*stk.wgt.error[,1:ny]
    return(mwgt.obs)
}

# Obs.fec
#---------------
# Function to simulate the fecundity observation at age (yearly mean), taking into account the age misclassification
# - Fec at age observation error
#Input
# biol       : an an object of class FLBiol
# ages.error : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# fec.error  : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# yr         : integer, the year the stock is observed from
#Details
#! Observation error of proportions are not additive nor multiplicative. Thus to introduce observation error in proportion
#!mature it is necessary to replace the real values in totto with the observed values.
#!In creating the replacement observed values the mean will be given by the real values and the shape2 parameter of the Beta distribution.
#!bot.age and top.age allows to constrain the uncertainty to a range of intermediate ages, taking advantage of the fact that very young fish 
#!(i.e. age 0) will not be reproducing no matter what, and obviously adult fish will be sexually mature (ceteris paribus).
#!If bot.age and top.age are equal then maturity will go from 0 at age (bot.age - 1) to 1 at age (top.age)
#! ACLARAR CON DORLETA
Obs.fec  <- function(biol, ages.error, fec.error, yr){
    
    ny   <- yr -1
    
    fec.real  <- seasonMeans(unitMeans(biol@fec))[,1:ny]
    fec.obs  <- fec.real
    it <- dim(fec.real)[6]
    
    for(i in 1:it){
        for(y in 1:ny){
            fec.obs[,y,,,,i] <- fec.real[,y,,,,i]%*%ages.error[,,y,i] 
        }
    }
    fec.obs  <- fec.obs*fec.error[,1:ny]
    
    return(fec.obs)
}

## OPERATING ON OBJECTS OF CLASS FLFleetExt
#############################################

# Obs.land.nage
#---------------
# Function to simulate the landings observation in numbers at age (year total), taking into account the age misclassification and the total error in landings
# - Landings in numbers at age observation error
#Input
# fleets          : an object of class FLFleetsExt
# ages.error      : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# land.nage.error : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# land.wgt.obs    : observed landings weight at age
# yr              : integer, the year the stock is observed from
# stknm           : character, name of the stock
#Details
# If land.nage.error = 1 numbers at age estimation error in landings is suppressed.
Obs.land.nage <- function(fleets, ages.error, land.nage.error, land.wgt.obs, yr, stknm){

    ny     <- yr -1

    la.num <- unitSums(seasonSums(landStock(fleets,stknm)))[,1:ny] #sum thru units of numbers of all season

    na     <- dim(la.num)[1]
    it     <- dim(la.num)[6]     
       
    ola.num <- array(0,c(na,ny,1,1,1,it))
    mwgt.real <- seasonMeans(unitMeans(wtalStock(fleets,stknm)))[,1:ny]

    for(j in 1:ny){
       for(n in 1:it){
           ola.num[,j,,,,n]  <- la.num[,j,,,,n]%*%ages.error[,,j,n]
        }
    }
     
    # Correct the landings at age so that _new_ total landings are equal to 
    # _real_ total landings (corrected by dga)
    Lnew  <- apply(ola.num*land.wgt.obs, c(2:6), sum)  # [ny,1,1,1,it] (the factors of the * are arrays)
    Lreal <- apply(la.num*mwgt.real, c(2:6), sum)  # [1,ny,1,1,1,it]  (the factors of the * are FLQ-s)
    Lrat <- Lreal/array(Lnew, dim = c(1, dim(Lnew)))
    
    Lrat[is.na(Lrat)] <- 0
    
    ola.num <- sweep(ola.num, 2:6, Lrat, "*")
    
    # Apply the error in the observation of the landings due to misreporting.
    ola.num          <- ola.num*land.nage.error[,1:ny]
    
    return(ola.num)
}

# Obs.land.wgt
#---------------
# Function to simulate the observation of the mean weight at age in the landings (seasonal mean), taking into account the age misclassification
# - Weight at age in landings observation error
#Input
# fleets         : an object of class FLFleetsExt
# ages.error     : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# land.wgt.error : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# land.nage.obs  : observed landings numbers at age (not currently used)
# yr             : integer, the year the stock is observed from
# stknm          : character, name of the stock
#Details
# If land.wgt.error = 1 weigth at age estimation error in landings is suppressed.
Obs.land.wgt <- function(fleets, ages.error, land.wgt.error, yr, stknm){

     ny   <- yr -1
     laa <- unitSums(seasonSums(landStock(fleets,stknm)))[,1:ny]
     
     it <- dim(fleets[[1]]@effort)[6]
     
     mwgt.real <- seasonMeans(unitMeans(wtalStock(fleets,stknm)))[,1:ny]
     mwgt.obs  <- mwgt.real
     
     for(i in 1:it){
        for(y in 1:ny){
       # [rr]     mwgt.obs[,y,,,,i] <- ((mwgt.real[,y,,,,i,drop=T]*laa[,y,,,,i,drop=T])%*%ages.error[,,y,i])/land.nage.obs[,y,,,,i, drop=T]
            mwgt.obs[,y,,,,i] <- mwgt.real[,y,,,,i,drop=T]%*%ages.error[,,y,i]

        }
     }
    
     mwgt.obs <- mwgt.obs*land.wgt.error[,1:ny]
     return(mwgt.obs)
}

# Obs.disc.nage
#---------------
# Function to simulate the discards observation in numbers at age (year total), taking into account the age misclassification and the total error in discards
# - Discards in numbers at age observation error
#Input
# fleets          : an object of class FLFleetsExt
# ages.error      : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# disc.nage.error : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# disc.wgt.obs    : observed discards weight at age
# yr              : integer, the year the stock is observed from
# stknm           : character, name of the stock
#Details
# If disc.nage.error = 1 numbers at age estimation error in discards is suppressed.
Obs.disc.nage <- function(fleets, ages.error, disc.nage.error, disc.wgt.obs, yr, stknm){
  
    ny     <- yr -1 
    da.num <- unitSums(seasonSums(discnStock(fleets,stknm)))[,1:ny] #sum thru units of numbers of all season
    it     <- dim(da.num)[6]
    na     <- dim(da.num)[1]
    
    oda.num <- da.num; oda.num[] <- 0
    
    mwgt.real <- seasonMeans(unitMeans(wtadStock(fleets,stknm)))[,1:ny]

    
    for(j in 1:ny){
        for(n in 1:it){
            oda.num[,j,,,,n] <- da.num[,j,,,,n]%*%ages.error[,,j,n]
        }
    }
    
         
    # Correct the discards at age so that _new_ total discards are equal to 
    # _real_ total discards (corrected by dga)
    Dnew  <- apply(oda.num*disc.wgt.obs, c(2:6), sum)  # [1,ny,1,1,1,it]
    Dreal <- apply(da.num*mwgt.real, c(2:6), sum)  # [1,ny,1,1,1,it] 
    Drat <- Dreal/Dnew

    Drat[is.na(Drat),] <- 1  # Drat == NA => is because Dnew = Dreal == 0.
    
    oda.num <- sweep(oda.num, 2:6, Drat, "*")
    
    # Apply the error in the observation of the discards due to misreporting.
    oda.num          <- oda.num*disc.nage.error[,1:ny]
    
    return(oda.num)
}

# Obs.disc.wgt
#---------------
# Function to simulate the observation of the mean weight at age in the discards (seasonal mean), taking into account the age misclassification
# - Weight at age in discards observation error
#Input
# fleets         : an object of class FLFleetsExt
# ages.error     : an array of probabilities of age assignment of dim = c(na, na, ny, iter)
# disc.wgt.error : an array of random multiplicative deviates of dim c(na,ny,1,1,1,it)
# disc.nage.obs  : observed landings numbers at age (not currently used)
# yr             : integer, the year the stock is observed from
# stknm          : character, name of the stock
#Details
# If disc.wgt.error = 1 mean weight at age estimation error in discards is suppressed.
Obs.disc.wgt <- function(fleets, ages.error, disc.wgt.error,  yr, stknm){
    
    ny   <- yr -1
    daa <- unitSums(seasonSums(discStock(fleets,stknm)))[,1:ny]
    
    it <- dim(fleets[[1]]@effort)[6]
    
    mwgt.real <- seasonMeans(unitMeans(wtadStock(fleets,stknm)))[,1:ny]
    mwgt.obs  <- mwgt.real
    
    for(i in 1:it){
       for(y in 1:ny){
      # [rr]     mwgt.obs[,y,,,,i] <- ((mwgt.real[,y,,,,i]*daa[,y,,,,i])%*%ages.error[,,y,i])/disc.nage.obs[,y,,,,i,drop=T]
            mwgt.obs[,y,,,,i] <- mwgt.real[,y,,,,i]%*%ages.error[,,y,i]
       }
    }
    
    mwgt.obs <- mwgt.obs*disc.wgt.error[,1:ny]
    
    return(mwgt.obs)
}

#Harvest at age from solving Baranov eq. with uniroot
Obs.FSolver <- function(Fay,May,Nay,Cay)
       {
       Cay-(Fay/(Fay+May))*(1-exp(-Fay-May))*Nay
       }


################################################################################
############ GLOBAL/AGE-STRUCTURED POP OBSERVED BY GLOBAL BIOMASS ##############
################################################################################

## OPERATING ON OBJECTS OF CLASS FLBiols
##########################################

# Obs.stk.bio
#---------------
# Function to simulate the observation of the total biomass in the population at the beginning of the year
# - Total biomass observation error
#Input
# biol          : an object of class FLBiols
# stk.bio.error : an array of multiplicative errors (total abundance estimation error), dim = c(1,ny,1,1,1,it)
#                 for systematic and/or random error in the observation of total biomass
# yr            : integer, the year of assessment/management
Obs.stk.bio <- function(biol, stk.bio.error, yr){
     btot     <- unitSums(biol@n[,,,1,,]*biol@wt[,,,1,,])
     ny       <- yr - 1 
     btot     <- btot*stk.bio.error[,1:ny]
     return(btot)
}

## OPERATING ON OBJECTS OF CLASS FLFleetExt
#############################################

# Obs.land.bio
#---------------
# Function to simulate the observation of the total landings with total error
# - Total landings observation error
#Input
# fleets         : an object of class FLFleetsExt
# land.bio.error : an array of multiplicative errors (total landings estimation error), dim = c(1,ny,1,1,1,it)
# yr             : integer, the year the stock is observed from
# stknm          : character, name of the stock
#Details
# If land.bio.error = 1 the observation of the total landings is perfect.
Obs.land.bio <- function(fleets, land.bio.error,  yr, stknm){                                                                       
    ny        <- yr - 1
    tland     <- unitSums(seasonSums(tlandStock(fleets, stknm)))[,1:ny]
    tland     <- tland*land.bio.error[,1:ny]

    return(tland)
}

# Obs.disc.bio
#---------------
# Function to simulate the observation of the total discards with total error
# - Total discards observation error
#Input
# fleets         : an object of class FLFleetsExt
# disc.bio.error : an array of multiplicative errors (total discards estimation error), dim = c(1,ny,1,1,1,it)
# yr             : integer, the year the stock is observed from
# stknm          : character, name of the stock
#Details
# If disc.bio.error = 1 the observation of the total discards is perfect.
Obs.disc.bio <- function(fleets, disc.bio.error, yr, stknm){
     ny        <- yr -1
     tdisc     <- unitSums(seasonSums(tdiscStock(fleets, stknm)))[,1:ny]
     tdisc     <- tdisc*disc.bio.error[,1:ny]

     return(tdisc)
}

#END