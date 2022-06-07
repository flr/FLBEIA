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
 #' @param fleets.ctrl Fleets' control file containing the catch production function for each fleet and stock.
 #' @param mean.yrs A character vector with the name of the years used to calculate mean selectivity.
 #' @param sim.yrs A character vector with the name of the years in the projection period.
 #' 
 #' @return A FLFleetsExt object. 
 #' 
 
 calculate.q.sel.flrObjs <- function(biols, fleets, BDs, fleets.ctrl, mean.yrs, sim.yrs){
   
    for(st in names(biols)){
    
      na <- dim(biols[[st]]@n)[1]
      
    # For age structured models calculate always Fa, to save on computations, for years 1:simyrs[]-1
      yrs <- dimnames(biols[[st]]@n)[[2]]

      if(na > 1) Fa <- Fa_cond_Baranov(biols, fleets, years = yrs, stk = st)
      Ct <- catchStock(fleets,st)
    
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
          
          cobj <- fleets[[fl]]@metiers[[mt]]@catches[[st]]
          
        #  if(fl == 'GN7_SP' & st == 'HKE') browser()
          
            if(!(st %in% catchNames(fleets[[fl]]@metiers[[mt]]))) next  
          
            catchProd <- fleets.ctrl[[fl]][[st]][['catch.model']] 
              
            C     <- (cobj@discards.n + cobj@landings.n)
            alpha <- cobj@alpha
            beta  <- cobj@beta
            E     <- fleets[[fl]]@effort*fleets[[fl]]@metiers[[mt]]@effshare
            
            if(na == 1 ) C <- cobj@discards.n*cobj@discards.wt + 
              cobj@landings.n*cobj@landings.wt
            
            
            # Cobb-Douglas q
            if(substr(catchProd,1,11) == 'CobbDouglas') cobj@catch.q <- C/((E%^%alpha)*(B%^%beta))
            # Baranov q: First calculate Fa for the whole fishery and then using partial F calculate q = Fap/E^alpha.
            if(catchProd == 'Baranov'){ 
              Fpa <- (C/Ct)*Fa
              
              cobj@catch.q[] <- Fpa/(E%^%alpha)
            }  
           
            cobj@landings.sel <- cobj@landings.n/(cobj@landings.n +
                                                                                                                             fleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n)
            cobj@discards.sel <- 1 - cobj@landings.sel
            
            # Fill in the values in the projection.
            cobj@catch.q[,ac(sim.yrs)] <- yearMeans(cobj@catch.q[,ac(mean.yrs)])
            cobj@landings.sel[,ac(sim.yrs)] <- apply(cobj@landings.sel[,ac(mean.yrs)],1,mean,na.rm = TRUE)
            cobj@discards.sel[,ac(sim.yrs)] <- 1-yearMeans(cobj@landings.sel[,ac(mean.yrs)])

            # If there are NA-s replace them by 0 in the case of catch.q & discards.sel and by 1 in the case of landings
            cobj@catch.q[is.na(cobj@catch.q)] <- 0
            cobj@landings.sel[is.na(cobj@landings.sel)] <- 1
            cobj@discards.sel[is.na(cobj@discards.sel)] <- 0
            
            mt_idx<-which(names(fleets[[fl]]@metiers)==mt)
            fleets[[fl]]<-fill_flcatches(fl=fleets[[fl]],cobj=cobj,st=st,mt_idx=mt_idx)
            rm(cobj)
            
            }
     }
   }
   return(fleets)
}
   
   
# Auxiliar function to calculate F-at-age.
 
 Fa_cond_Baranov <-    function(biols, fleets, stk, years = dimnames(biols[[1]]@n)$year){
     stknms <- names(biols)
     
     it     <- dim(biols[[1]]@n)[6]
     ny     <- length(years)
     ns     <- dim(biols[[1]]@n)[4]
     yrnms  <- years
     ssnms <- dimnames(biols[[1]]@n)[[4]]
     
     # harvest: * if age structured calculate it from 'n'.
       #          * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
       na <- dim(biols[[stk]]@n)[1]
       nu <- dim(biols[[stk]]@n)[3]
       
       res <- unclass(biols[[stk]]@n)
       res[] <- NA
       
       Dnms <- dimnames(biols[[stk]]@n)
       
       n.  <- unclass(biols[[stk]]@n[,years,,,drop=F])
       m.  <- unclass(biols[[stk]]@m[,years,,,drop=F])
       c.  <- unclass(catchStock(fleets, stk)[,years,,,drop = F])
       
       fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
       
       for(ss in ssnms){
         for(y in yrnms){
           for(a in 1:na){
             for(u in 1:nu){
               for(i in 1:it){
                 if(is.na(n.[a,y,u,ss,1,i])){ 
                   res[a,y,u,ss,1,i] <- NA
                 }
                 else{
                   if(n.[a,y,u,ss,1,i] == 0){ res[a,y,u,ss,1,i] <- 0; next}
                   # if n. < c. we take the Fa in the previous 'correct  age because otherwise it could generate problems
                   if(n.[a,y,u,ss,1,i] < c.[a,y,u,ss,1,i]){
                     if(a == 1) res[a,y,u,ss,1,i] <- Inf
                     else res[a,y,u,ss,1,i] <- res[a-1,y,u,ss,1,i]}
                   else{
                     xx <- try(uniroot(fobj, lower = 0, upper = 1e6, n = n.[a,y,u,ss,1,i], m=m.[a,y,u,ss,1,i], c = c.[a,y,u,ss,1,i])$root, silent = TRUE)
                     res[a,y,u,ss,1,i] <- ifelse(class(xx) == 'try-error', NA, xx)
                   }
                 }     
               }}}}}
     
     return(res)
   }
   
  
