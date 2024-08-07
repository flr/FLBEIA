#-------------------------------------------------------------------------------
#                   POPULATION GROWTH FUNCTIONS             
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMARK: '...' in the arguments of the functions are neccesary in order to be
#   generalistic inside 'biols.om' function ('eval(call(...)').
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#------------------------------------------------
# ** ASPG(biol, SR, fleets, biol.control) ~  Age Structured Population Growth.
# Projects in one season an age structured population using
#    << N[s] = (N[s-1]*exp(-M[s-1]/2) - Catch[s-1])*exp(-M[s-1]/2)  >>
# The assumption is that natural mortality is constant and continious 
# in each season but the catch (fishing mortality) occurs instantaneously 
# in the middle of 
# the season. 
#------------------------------------------------
#
#------------------------------------------------
# # ** BDPG(biol, BD, fleets, biol.control, ...) ~  Biomass Dynamic Population Growth.
# Projects in one season a biomass dynamic population using
#    << B[s] = B[s-1] - C[s-1] + Production[s-1]  >>
# Production is determined using BDsim. 
#------------------------------------------------
#
# Dorleta Garcia
# Created: 21/10/2010 10:08:16
# Changed: 27/10/2010 08:34:17
#-------------------------------------------------------------------------------

fixedPopulation <- function(biols, SRs, fleets, year, season, stknm, ...)  return(list(biol = biols[[stknm]]))

#-------------------------------------------------------------------------------
# ASPG(biol, SR, fleets, biol.control)
# - OUTPUT: list(biol = biol, SR = SR) - Updated FLBiol and FLSRsim objects.
#-------------------------------------------------------------------------------

ASPG <- function(biols, SRs, fleets, year, season, stknm, ...){

    cat('-----------------ASPG-----------\n')

    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    biol <- biols[[stknm]]     
    SR   <- SRs[[stknm]]
        
    dimnms <- dimnames(biol@n)
    
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
    
    na <- dim(biol@n)[1]
    ns <- dim(biol@n)[4]
    stock <- biol@name
    
    # IF season = 1 THEN age groups move to the next. 
    if(ss == 1){
        # total catch in year [y-1] season [ns].
        catch.n <- catchStock(fleets,stock)[,yr-1,,ns]
        # middle ages
        biol@n[-c(1,na),yr,,ss] <- (biol@n[-c(na-1,na),yr-1,,ns]*exp(-biol@m[-c(na-1,na),yr-1,,ns]/2) - catch.n[-c(na-1,na),])*
                                            exp(-biol@m[-c(na-1,na),yr-1,,ns]/2) 
        # plusgroup
        biol@n[na,yr,,ss]       <- (biol@n[na-1,yr-1,,ns]*exp(-biol@m[na-1,yr-1,,ns]/2) - catch.n[na-1,])*exp(-biol@m[na-1,yr-1,,ns]/2) + 
                                   (biol@n[na,yr-1,,ns]*exp(-biol@m[na,yr-1,,ns]/2) - catch.n[na,])*exp(-biol@m[na,yr-1,,ns]/2)

    }
    else{
        # total catch in year [yr] season [ss-1].
        catch.n <- catchStock(fleets,stock)[,yr,,ss-1]
        # middle ages      # for unit == ss  and age = 1, it will be equal NA but be updated after with SRsim.
        biol@n[,yr,,ss] <- (biol@n[,yr,,ss-1]*exp(-biol@m[,yr,,ss-1]/2) - catch.n)*exp(-biol@m[,yr,,ss-1]/2) 
    }

    # Update SSB.
    SR@ssb[,yr,,ss] <- unitSums(quantSums(n(biol) * wt(biol) * fec(biol)*mat(biol) * 
                                            exp(-biol@m*spwn(biol)), na.rm=TRUE))[,yr,,ss]

    # RECRUITMENT
    if(dim(biol@n)[3] > 1 & dim(biol@n)[3] == dim(biol@n)[4]){
        # 'number_of_units > 1' => Recruitment occurs in all seasons, the 'unit' correspond with the recruitment 'season'.
        SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
        biol@n[1,yr,ss,ss] <- SR@rec[,yr,,ss]
        biol@n[1,yr,-(1:ss),ss] <- 0  # The recruitment is 0 in [units > ss].
    }
    else{  # dim(biol@n)[3] = 1, The recruitment only occurs in 1 season every year. 
        if(SR@proportion[,yr,,ss,,1] == 1){ # If the recruitment season is 'ss' generate it otherwise
            SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
            biol@n[1,yr,,ss] <- SR@rec[,yr,,ss]
       } else if (ss==1) { # If the recruitment season is NOT the first one, set numbers at first age group to 0 in this season.
      biol@n[1,yr,,ss] <- 0
    } # If the recruitment season is NOT the first one, do nothing, the population in first age group is just the survivors of previous season.
    
    
    }
    
    if(any(biol@n[,yr,,ss]<0)){
        biol <- correct.biomass.ASPG(biol, yr, ss)
    }
    
    return(list(biol = biol, SR = SR))

} 


#-------------------------------------------------------------------------------
# BDPG(biol, BD, fleets, biol.control)
# - OUTPUT: list(biol = biol, BD = BD) - Upadated FLBiol and FLSRsim objects.
#-------------------------------------------------------------------------------

BDPG <- function(biols, BDs, fleets, year, season, stknm, ...){

    cat('-----------------BDPG-----------\n')

    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    biol <- biols[[stknm]]    
    BD   <- BDs[[stknm]]
    
    dimnms <- dimnames(biol@n)
    
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
    
    ns <- dim(biol@n)[4]
    stock <- biol@name
    
    # IF season = 1 
    if(ss == 1){
        yr0 <- yr - 1
        ss0 <- ns

    }
    else{
        yr0 <- yr
        ss0 <- ss - 1 
    }
    
    # total catch in year [yr0] season [ss0].
    BD@catch[,yr0,,ss0] <- catchStock(fleets,stock)[,yr0,,ss0]
        
    # Update FLBDsim object.
    BD <- BDsim(BD,yr,ss)
        
    # Update FLBiol Object
    biol@n[,yr,,ss] <- BD@biomass[,yr,,ss]/biol@wt[,yr,,ss]
    
    return(list(biol = biol, BD = BD))

} 

# correct.biomass.ASPG <- function(biol)
correct.biomass.ASPG <- function(biol, year, season){

    for(i in 1:dim(biol@n)[6]){
        N  <- c(biol@n[,year,, season,,i])
        wt <- c(biol@wt[,year,,season,,i])
    
        # identify the ages that are < 0, this should be only happen in the 
        # first year of simulation, when catch has not been calculated using CobbDoug.
        negs <- which(N<0)
        #biomass proportion in the positive ages, the negative biomass is discounted proportionally 
        # depending on the abundace of each age group. The calculation is done for all cohorts (units) at the same time.
        p <- N[-negs]*wt[-negs]/sum(N[-negs]*wt[-negs])
        # negative biomass
        NegB <- -sum(N[negs]*wt[negs])
        # New Biomass
        Nnew <- N                                       
        Nnew[negs]  <- 0
        Nnew[-negs] <- N[-negs] - (NegB*p/wt[-negs]) 
        
        biol@n[,year,, season,,i] <- matrix(Nnew, dim(biol@n)[1])
        }
        return(biol)
        
}
        


#-------------------------------------------------------------------------------
# The same as ASPG but the catch is produced using Baranov catch equation, i.e, 
#  the catch is done all along the year using an instantaneous constant 
#  fishing mortality rate. 
#-------------------------------------------------------------------------------

ASPG_Baranov <- function(biols, SRs, fleets, year, season, stknm, ...){
  
  cat('-----------------ASPG_Baranov-----------\n')
  
  if(length(year) > 1 | length(season) > 1)
    stop('Only one year and season is allowed' )
  
  biol <- biols[[stknm]]     
  SR   <- SRs[[stknm]]
  
  dimnms <- dimnames(biol@n)
  
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
  
  na <- dim(biol@n)[1]
  ns <- dim(biol@n)[4]
  ni <- dim(biol@n)[6]
  nu <- dim(biol@n)[3]
  stock <- biol@name
  

  # IF season = 1 THEN age groups move to the next. 
  if(ss == 1){
    # total catch in year [y-1] season [ns].
    catch.n <- catchStock(fleets,stock)[,yr-1,,ns]
    
    Ca <-     catch.n <- catchStock(fleets,stock)[,yr-1,,ns]
    
    findF <- function(Fa,Ca, Ma, Na){
      Ca. <- (Fa/(Fa+Ma))*(1-exp(-Ma-Fa))*Na
      res <- Ca.-Ca
      return(res)
    }
    
    
    catch.n <- catchStock(fleets,stock)[,yr-1,,ns]
    Ma <- unname(unclass(biol@m[,yr-1,,ns]))
    Na <- unname(unclass(biol@n[,yr-1,,ns]))
    Ca <- unname(unclass(catch.n))
    fa <- Ma # the same dimensions as Ma
    fa[] <- NA
    
    for(u in 1:nu){
      
    loop.uniroot <- function(i) {
      if(Ca[a,,u,,,i] > Na[a,,u,,,i]) return(10)
      if(Ca[a,,u,,,i] == 0) return(0) # Ca=0 --> Fa=0 (to avoid Inf when Ma=0)
      return(uniroot(findF,interval=c(0,2),Ca=Ca[a,,u,,,i],Ma=Ma[a,,u,,,i], Na=Na[a,,u,,,i], tol = 1e-12,extendInt = "yes")$root)
    }

    for (a in 1:na) fa[a,,u,,,] <- vapply(1:ni, loop.uniroot, numeric(1))
 
    za <- Ma+fa
    
    # middle ages
    biol@n[-c(1,na),yr,u,ss] <- biol@n[-c(na-1,na),yr-1,u,ns]*exp(-za[-c(na-1,na),,u,,,,drop=F])
    # plusgroup
    biol@n[na,yr,u,ss]       <- biol@n[na-1,yr-1,u,ns]*exp(-za[na-1,,u,,,,drop=F])+biol@n[na,yr-1,u,ns]*exp(-za[na,,u,,,,drop=F])
    # 
    }
  }
  else{
    # total catch in year [yr] season [ss-1].
    
    findF <- function(Fa,Ca, Ma, Na){
      Ca. <- (Fa/(Fa+Ma))*(1-exp(-Ma-Fa))*Na
      res <- sum(Ca.)-sum(Ca)
      return(res)
    }
    
    catch.n <- catchStock(fleets,stock)[,yr,,ss-1]
    Ma <- unname(unclass(biol@m[,yr,,ss-1]))
    Na <- unname(unclass(biol@n[,yr,,ss-1]))
    Ca <- unname(unclass(catch.n))
    fa <- Ma # the same dimensions as Ma
    fa[] <- NA
    
    for(u in 1:nu){
      
    loop.uniroot <- function(i) {
      if(Ca[a,,u,,,i] > Na[a,,u,,,i]) return(10)
      if(Ca[a,,u,,,i] == 0) return(0) # Ca=0 --> Fa=0 (to avoid Inf when Ma=0)
      uniroot(findF,interval=c(0,2),Ca=Ca[a,,u,,,i],Ma=Ma[a,,u,,,i], Na=Na[a,,u,,,i], tol = 1e-12,extendInt = "yes")$root
    }
    
    for (a in 1:na) fa[a,,u,,,] <- vapply(1:ni, loop.uniroot, numeric(1))
    
    za <- Ma+fa
    
    # middle ages      # for unit == ss  and age = 1, it will be equal NA but be updated after with SRsim.
    biol@n[,yr,u,ss] <- biol@n[,yr,u,ss-1]*exp(-za[,,u,,,,drop=F])
  }}
  
  # Update SSB.
  SR@ssb[,yr,,ss] <- unitSums(quantSums(n(biol) * wt(biol) * fec(biol)*mat(biol) * 
                                          exp(-biol@m*spwn(biol)), na.rm=TRUE))[,yr,,ss]
  
  # RECRUITMENT
  if(dim(biol@n)[3] > 1 & dim(biol@n)[3] == dim(biol@n)[4]){
    # 'number_of_units > 1' => Recruitment occurs in all seasons, the 'unit' correspond with the recruitment 'season'.
    SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
    biol@n[1,yr,ss,ss] <- SR@rec[,yr,,ss]
    biol@n[1,yr,-(1:ss),ss] <- 0  # The recruitment is 0 in [units > ss].
  }
  else{  # dim(biol@n)[3] = 1, The recruitment only occurs in 1 season every year. 
    if(SR@proportion[,yr,,ss,,1] == 1){ # If the recruitment season is 'ss' generate it otherwise
      SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
      biol@n[1,yr,,ss] <- SR@rec[,yr,,ss]
    } else if (ss==1) { # If the recruitment season is NOT the first one, set numbers at first age group to 0 in this season.
      biol@n[1,yr,,ss] <- 0
    } # If the recruitment season is NOT the first one, do nothing, the population in first age group is just the survivors of previous season.
    
    
  }
  
  if(any(biol@n[,yr,,ss]<0)){
    biol <- correct.biomass.ASPG(biol, yr, ss)
  }
  
  return(list(biol = biol, SR = SR))
  
} 

    
#-------------------------------------------------------------------------------
# ASPG_DDW(biol, SR, fleets, biol.control)
# - OUTPUT: list(biol = biol, SR = SR, covars = covars) - Updated FLBiol and FLSRsim objects.
#-------------------------------------------------------------------------------

# Age structured population growth with densodependence

ASPG_DDW <- function(biols, SRs, fleets, stknm, year, season, ctrl, covars, ...){
  
  cat('-----------------ASPG_DDW-----------\n')
  
  # check if DDW covars available
  if (!"DDW" %in% names(covars) | !stknm %in% names(covars[["DDW"]]))
    covars[["DDW"]][[stknm]] <- biols[[stknm]]@n * 0 + 1
  
  # estimate the DD weights
  ddw.model <- ctrl[[stknm]][['ddw.model']]
  ddw.ctrl  <- ctrl[[stknm]][['ddw.ctrl']]
  
  biol <- biols[[stknm]]
  
  if (ddw.model == "ddwAgeCa") # for in-year ssb or total biomass calculation (nage required)
    biol <- ASPG(biols, SRs, fleets, year, season, stknm, biols.ctrl,...)$biol
 
  wts <- eval(call(ddw.model, biol = biol, stknm = stknm, year = year, season = season, 
                  ctrl = ddw.ctrl, covars = covars))
  
  # save
  biols[[stknm]]@wt[,year,,season,]        <- wts$wt
  covars[["DDW"]][[stknm]][,year,,season,] <- wts$wt.chg
  
  # use normal ASPG to project the population
  
  res <- ASPG(biols, SRs, fleets, year, season, stknm, biols.ctrl,...)
  
  res$covars <- covars
  
  return(res)
}
