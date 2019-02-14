###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:       calculate.CDparam
# NOTE #1:     calculates Cobb Douglass parameters: 
#              returns a list with alpha, beta and catch.q. Alpha and beta are 
#              FLQuants equal one and 
#              catch.q is: catch/(effshare*abundance)
###############################################################################
#-------------------------------------------------------------------------
#  inputs:
#   flq:  FLquant with the dimensions of the stock and age='all'
#   flqa: FLquant with the dimensions of the stock and age from min to max.
#   age.min: Minimum age of the stock
#   age.max: Maximum age of the stock
#   landings.n: Landings at age (FLQuant	[na,ny(hist),nu(stock),ns,1/ni)])
#   discards.n: Discards at age (FLQuant  [na,ny(hist),nu(stock),ns,1/ni)])
#   effort:	  fleet's effort (FLQuant	[1,ny(hist),nu(stock),ns,1/ni)])
#   effshare: effort share per metier (FlQuant	[1,ny(hist),nu(stock),ns,1/ni)])
#   stk.n:    Stock abundance at age
# Optional:
#   largs: Other arguments. For stocks in biomass largs$stk.gB: surplus production
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        Check inputs
#   Section 2:        Set dimensions
#   Section 3:        Calculate catch.q
#   Section 4:        Set alpha, beta
#   Section 5:        Create a list with the 3 parameters
#-------------------------------------------------------------------------------

 calculate.CDparam <- function(stk.n,landings.n,discards.n,
                                  effort,effshare,age.min,age.max,flqa,flq,largs) { 

  #==============================================================================
  # Section 1:        Check inputs
  #==============================================================================
  
  if(all(is.na(effort)) || all(effort==0)) {
    stop('Cobb Douglas parameters: Effort data missing')}

  if(all(is.na(effshare)) || all(effshare==0)) {
    stop('Cobb Douglas parameters: Effort share data missing')}
  
  # if(any(is.na(age.min)|| is.na(age.max))) {
  #   stop('Cobb Douglas parameters: Cobb Douglas parameters: Wrong min or max age')}
  
  if(all(is.na(landings.n)) || all(is.na(discards.n))) {
    stop('Cobb Douglas parameters: Na-s in landings.n and discards.n')}
   if(all(is.na(discards.n))) {
     warning('warning: all NA-s in discards.n')}
  
  #==============================================================================
  # Section 2:        Set dimensions
  #==============================================================================
  
  catch.q <- flqa
  alpha   <- flqa
  beta    <- flqa

  hist.yrs     <- dimnames(effort)$year
  met.effort   <- flq[,hist.yrs]
  met.effort[,hist.yrs] <- effort*effshare        
    
  #==============================================================================
  # Section 3:        Calculate catch.q
  #==============================================================================
  if(all(is.na(c(age.max,age.min)))){
    stk.gB <- largs$stk.gB
    catch.q[,hist.yrs]   <- ((landings.n + discards.n)[,hist.yrs])/(met.effort*(stk.n+stk.gB[,hist.yrs]))
  }else{
    if(length(age.min:age.max)==1 ){
     stk.gB <- largs$stk.gB
     catch.q[,hist.yrs]   <- ((landings.n + discards.n)[,hist.yrs])/(met.effort*(stk.n+stk.gB[,hist.yrs]))
    }else{
    for (aa in 1:length(age.min:age.max)){ 
      stk.m <- largs$stk.m
      catch.q[aa,hist.yrs]   <- ((landings.n + discards.n)[aa,hist.yrs])/(met.effort*stk.n[aa,]*exp(-stk.m[aa,]/2))}
    }}       
 
   catch.q[is.infinite(catch.q)] <- 0
   catch.q[is.na(catch.q)] <- 0
  #==============================================================================
  # Section 4:         Set alpha, beta
  #==============================================================================
  
  alpha[]   <- 1
  beta[]    <- 1

  #==============================================================================
  # Section 5:         Create a list with the 3 parameters
  #==============================================================================
  
  CD_param <- vector('list',3)
  names(CD_param) <- c('alpha','beta','catch.q')
  CD_param[['alpha']]    <- alpha
  CD_param[['beta']]     <- beta
  CD_param[['catch.q']]  <- catch.q
  
  return(CD_param)
}
  

 
 #==============================================================================
 # Section 6:       
 #==============================================================================
 
 
 #' Calculates selectivity-related parameters
 #' 
 #' Given stock abundances and catches (i.e. landings and discards), this function estimates values for
 #' \code{fleets[[fl]]@@metiers[[mt]]@@catches[[st]]@@catch.q}, 
 #' \code{fleets[[fl]]@@metiers[[mt]]@@catches[[st]]@@landings.sel}, and 
 #' \code{fleets[[fl]]@@metiers[[mt]]@@catches[[st]]@@discards.sel} 
 #' for all years except the simulation years. For the simulation years, 
 #' the parameter values are set as the mean of the parameters along \code{mean.yrs}. 
 #'
 #' @param biols An FLBiols object.
 #' @param fleets An FLFleetsExt object. An extended version of the FLFleet object defined in FLCore. 
 #' @param BDs A list of FLSRsim objects. One per biomass dynamic stock in biols object.
 #' @param mean.yrs A character vector with the name of the years used to calculate mean selectivity.
 #' @param sim.yrs A character vector with the name of the years in the projection period.
 #' 
 #' @return A FLFleetsExt object. 
 #' 
 
 calculate.q.sel.flrObjs <- function(biols, fleets, BDs, mean.yrs, sim.yrs){
   
    for(st in names(biols)){
    
      na <- dim(biols[[st]]@n)[1]
    
     if(na != 1){  # 'Biomass' in numbers because the catch is in numbers, in the middle of the season.
        B <- biols[[st]]@n*exp(-biols[[st]]@m/2)#*biols[[st]]@wt  
     }else{ # 'Biomass' in weight because the growth is in weight => later we use the catch in weight.
        
        if(is.null(BDs[[st]])) gB <- 0
        else gB <- BDs[[st]]@gB
        
        B <- biols[[st]]@n*biols[[st]]@wt + gB
     }
     
     for(fl in names(fleets)){
        for(mt in names(fleets[[fl]]@metiers)){
            
          cat(fl, ' - ', mt, ' - ', st, '\n')
            if(!(st %in% catchNames(fleets[[fl]]@metiers[[mt]]))) next  
          
            C     <- (fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n + fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n)
            alpha <- fleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha
            beta  <- fleets[[fl]]@metiers[[mt]]@catches[[st]]@beta
            E     <- fleets[[fl]]@effort*fleets[[fl]]@metiers[[mt]]@effshare
            
            if(na == 1 ) C <- fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n*fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt + 
                              fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n*fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt
            
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q <- C/((E%^%alpha)*(B%^%beta))
            
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel <- fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n/(fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n +
                                                                                                                             fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n)
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel <- 1 - fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel
            
            # Fill in the values in the projection.
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,ac(sim.yrs)] <- yearMeans(fleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,ac(mean.yrs)])
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,ac(sim.yrs)] <- apply(fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,ac(mean.yrs)],1,mean,na.rm = TRUE)
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel[,ac(sim.yrs)] <- 1-yearMeans(fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,ac(mean.yrs)])

            # If there are NA-s replace them by 0 in the case of catch.q & discards.sel and by 1 in the case of landings
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[is.na(fleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q)] <- 0
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[is.na(fleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel)] <- 1
            fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel[is.na(fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel)] <- 0
            }
     }
   }
   return(fleets)
}
   
   
   
   
  
