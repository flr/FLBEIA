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
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        Check inputs
#   Section 2:        Set dimensions
#   Section 3:        Calculate catch.q
#   Section 4:        Set alpha, beta
#   Section 5:        Create a list with the 3 parameters
#-------------------------------------------------------------------------------

 calculate.CDparam <- function(stk.n,landings.n,discards.n,
                                  effort,effshare,age.min,age.max,flqa,flq) { 

  #==============================================================================
  # Section 1:        Check inputs
  #==============================================================================
  
  if(all(is.na(effort)) || all(effort==0)) {
    stop('Cobb Douglas parameters: Effort data missing')}

  if(all(is.na(effshare)) || all(effshare==0)) {
    stop('Cobb Douglas parameters: Effort share data missing')}
  
  if(any(is.na(age.min)|| is.na(age.max))) {
    stop('Cobb Douglas parameters: Cobb Douglas parameters: Wrong min or max age')}
  
  if(all(is.na(landings.n))|| all(is.na(discards.n))) {
       stop('Cobb Douglas parameters: Na-s in landings.n and discards.n')}
  
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
  
  for (aa in 1:length(age.min:age.max)){  
    catch.q[aa,hist.yrs]   <- ((landings.n + discards.n)[aa,hist.yrs])/(met.effort*stk.n[aa,])}
        
  catch.q[is.infinite(catch.q)] <- 0
  
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
  
  
  
