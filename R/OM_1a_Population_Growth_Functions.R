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


#-------------------------------------------------------------------------------
# ASPG(biol, SR, fleets, biol.control)
# - OUTPUT: list(biol = biol, SR = SR) - Upadated FLBiol and FLSRsim objects.
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
    SR@ssb[,yr,,ss] <- unitSums(quantSums(n(biol) * wt(biol) * fec(biol) * exp(-spwn(biol) * m(biol)), na.rm=TRUE))[,yr,,ss]
        
    if(dim(biol@n)[3] > 1){
        # 'number_of_units > 1' => Recruitment occurs in all seasons, the 'unit' correspond with the recruitment 'season'.
        SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
        biol@n[1,yr,ss,ss] <- SR@rec[,yr,,ss]
        biol@n[1,yr,-(1:ss),ss] <- 0  # The recruitment is 0 in [units > ss].
    }
    else{  # dim(biol@n)[3] = 1, The recruitment only occurs in 1 season every year. 
        if(SR@proportion[,yr,,ss,,1] == 1){ 
            SR <- SRsim(SR, year = yr, season = ss, iter = 'all') 
            biol@n[1,yr,,ss] <- SR@rec[,yr,,ss]
        }
    
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
