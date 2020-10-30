#-------------------------------------------------------------------------------
#              - Functions to run basic checks.    
#
# Created: 17/01/2011 14:47:19
# Author: Dorleta Garcia
# Changed: 2019-03-15 15:57:45 - Sonia Sanchez (included more functions)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  checkDims(object)    - check that all the objects have the same dimension
#    in year ans season and it (it can be 1 or it but common it). 
#   minyear, maxyear: Characters.
#   ns, it: numeric
#-------------------------------------------------------------------------------


checkDims <- function(object, minyear, maxyear, ns, it){
  
  for(k in 1:length(object)){
    zzz <- unlist(dims(object[[k]]))
    if(minyear != zzz['minyear']) stop('Wrong "minyear" in', names(object)[k], ' element within ', class(object), ' object')
    if(maxyear != zzz['maxyear']) stop('Wrong "maxyear" in', names(object)[k], ' element within ', class(object), ' object')
    if(ns != as.numeric(zzz['season'])) stop('Wrong number of seasons in', names(object)[k], ' element within ', class(object), ' object')
    if(it != as.numeric(zzz['iter'])) stop('Wrong number of iterations in', names(object)[k], ' element within ', class(object), ' object')
  }
  return(TRUE)
} 



#' @title Checking inputs of FLBEIA function
#' @description This functions checks wether some conditions are met for the inputs 
#'              in the different arguments of the FLBEIA function.
#'              There are also checking functions for specific inputs: 
#'              biols, fleets, SRs, advice, and obs.ctrl arguments.
#' 
#' Functions available:
#' \itemize{
#' 
#'    \item{\code{checkFLBEIAData}}{: for checking inputs for differents arguments of FLBEIA function 
#'                               (of class FLBiols).}
#'     
#'    \item{\code{checkBiols}}{: for checking biols argument of FLBEIA function 
#'                               (of class FLBiols).}
#'    \item{\code{checkFleets}}{: for checking fleets argument of FLBEIA function 
#'                               (of class FLFleetsExt).}
#'    \item{\code{checkSRs}}{: for checking SRs argument of FLBEIA function 
#'                               (list of FLSRsim objects).}
#'    \item{\code{checkAdvice}}{: for checking advice argument of FLBEIA function 
#'                               (of class list with two FLQuant elements, TAC and quota.share).}
#'    \item{\code{checkObsctrl}}{: for checking obs.ctrl argument of FLBEIA function 
#'                               (of class list).}
#' }
#' 
#' @param biols       An FLBiols object.
#' @param fleets      An FLFleetsExt object.
#' @param SRs         A list of FLSRsim objects.
#' @param BDs         A list of FLBDsim objects.
#' @param covars      A list of FLQuants.
#' @param indices     A list of FLIndices.
#' @param advice      A list with two FLQuant elements, TAC and quota.share.
#' @param main.ctrl   A list with the settings to control the main function.
#' @param biols.ctrl  A list with the settings to control the biological operating model for each stock.
#' @param fleets.ctrl A list with the settings to control the fleets operating model for each fleet.
#' @param covars.ctrl A list with the settings to control the covars operating model for each fleet.
#' @param obs.ctrl    A list with the settings to control the observation model for each stock.
#' @param assess.ctrl A list with the settings to control the specify the assessment model for each stock.
#' @param advice.ctrl A list with the settings to control the advice model for each stock.
#' @param object      An object of the appropiate class depending on the FLBEIA argument to be checked:
#'                    \itemize{
#'                      \item{\code{biols}}   {: A FLBiols object.}
#'                      \item{\code{fleets}}  {:A FLFleetsExt object.}
#'                      \item{\code{SRs}}     {: A list of FLSRsim objects.}
#'                      \item{\code{BDs}}     {: A list of FLBDsim objects.}
#'                      \item{\code{advice}}  {: A list with two FLQuant elements, TAC and quota.share..}
#'                      \item{\code{obs.ctrl}}{: A list.}
#'                    }
#'
#'
#' @param ctrl   An object of class list (optional argument), input for fleets.ctrl argument of FLBEIA function.
#'
#' @return An error message if there is any error detected.
#' 



#-------------------------------------------------------------------------------
#  checkBiols(object) - checking biols argument of FLBEIA function
#-------------------------------------------------------------------------------

checkFLBEIAData <- function(biols, SRs = NULL, BDs = NULL, fleets, covars = NULL, indices = NULL, advice = NULL, main.ctrl, biols.ctrl, 
                            fleets.ctrl, covars.ctrl, obs.ctrl, assess.ctrl,  advice.ctrl){
  
  sim.years <- as.numeric(main.ctrl$sim.years)
  
  
  # Check input objects
  
  # - biols
  
  checkBiols(lapply(biols,window,start=sim.years[1]-1,end=sim.years[2]))
  
  # - SRs
  
  checkSRs(lapply(SRs,window, start=sim.years[1]-1,end=sim.years[2]))
  
  # - BDs
  
  checkBDs(lapply(BDs,window, start=sim.years[1]-1,end=sim.years[2]))
  
  # - fleets
  
  checkFleets(window(fleets,start=sim.years[1]-1,end=sim.years[2]), ctrl = fleets.ctrl)
  
  # - advice
  
  advice$quota.share <- lapply(advice$quota.share, window, start=sim.years[1]-1, end=sim.years[2])
  checkAdvice(advice)
  
  # - obs.ctrl
  
  checkObsctrl(obs.ctrl)
  
  
  # Extract the common dimensions [year, season, it] from the 1st biol.
  
  ny <- dim(biols[[1]]@n)[2]
  ns <- dim(biols[[1]]@n)[4]
  it <- dim(biols[[1]]@n)[6]
  minyear <- ac(dims(biols[[1]])$minyear)
  maxyear <- ac(dims(biols[[1]])$maxyear)
  
  
  # Check that all FLQuants have the rigth [ny,ns,it] dimensions. 
  
  chckdim0 <- checkDims(biols,  minyear, maxyear, ns, it)
  chckdim1 <- checkDims(fleets, minyear, maxyear, ns, it)
  if(!is.null(covars)) chckdim2 <- checkDims(covars, minyear, maxyear, ns, it)
  
  
  return(TRUE)
  
}


# data(one)
# checkFLBEIAData( biols = oneBio, SRs = oneSR, BDs = NULL, fleets = oneFl, 
#                  covars = oneCv, indices = NULL, advice = oneAdv, 
#                  main.ctrl = oneMainC, biols.ctrl = oneBioC, fleets.ctrl = oneFlC, 
#                  covars.ctrl = oneCvC, obs.ctrl = oneObsC, assess.ctrl = oneAssC, advice.ctrl = oneAdvC)
# data(oneIt)
# checkFLBEIAData( biols = oneItBio, SRs = oneItSR, BDs = NULL, fleets = oneItFl, 
#                  covars = oneItCv, indices = NULL, advice = oneItAdv, 
#                  main.ctrl = oneItMainC, biols.ctrl = oneItBioC, fleets.ctrl = oneItFlC, 
#                  covars.ctrl = oneItCvC, obs.ctrl = oneItObsC, assess.ctrl = oneItAssC, advice.ctrl = oneItAdvC)
# data(multi)
# checkFLBEIAData( biols = multiBio, SRs = multiSR, BDs = multiBD, fleets = multiFl, 
#                  covars = multiCv, indices = NULL, advice = multiAdv, 
#                  main.ctrl = multiMainC, biols.ctrl = multiBioC, fleets.ctrl = multiFlC, 
#                  covars.ctrl = multiCvC, obs.ctrl = multiObsC, assess.ctrl = multiAssC, advice.ctrl = multiAdvC)


#-------------------------------------------------------------------------------
#  checkBiols(object) - checking biols argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkFLBEIAData
#' @aliases checkBiols
 
checkBiols <- function(object) {
  
  for (st in names(object)) {
    
    # mat: 0 < maturity at age < 1
    if ( any(mat(object[[st]]) < 0) | any(mat(object[[st]]) > 1) ) 
      stop("Check mat for stock '", st,"' (required 0 < mat < 1)." )
    
    # spwn: 0 < proportion of mortality before spawning < 1
    if ( any(mat(object[[st]]) < 0) | any(mat(object[[st]]) > 1) ) 
      stop("Check spwn for stock '", st,"' (required 0 < spwn < 1)." )
    
  }
  
  return(TRUE)
  
}


# data(one)
# checkBiols(oneBio)
# data(oneIt)
# checkBiols(oneItBio)
# data(multi)
# checkBiols(multiBio)
# 
# # Errors
# obj1 <- obj2 <- oneBio
# mat(obj1$stk1)[1,1,] <- -0.5 # mat < 0
# checkBiols(obj1)
# mat(obj2$stk1)[1,1,] <- 5
# checkBiols(obj2)



#-------------------------------------------------------------------------------
#  checkFleets(object) - checking fleets argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkFLBEIAData
#' @aliases checkFleets

checkFleets <- function(object, ctrl=NULL) {
  
  for (fl in names(object)) {
    
    # capacity: only NA if "fixedEffort"
    if ( any(is.na(object[[fl]]@capacity)) )
      if (!is.null(ctrl)) {
        if ( ctrl[[fl]]$effort.model != "fixedEffort" ) 
          stop("Selected effort model for fleet '", fl,"' ('", 
               ctrl[[fl]]$effort.model,"') requires not NA values for capacity.")
      } else 
        warning("Check capacity if fleets.ctrl[['",fl,"']]$effort.model != 'fixedEffort'")
    
    # effshare: sum(effort share by metier)==1, if effort>0
    if ( any(round(Reduce("+", lapply( object[[fl]]@metiers, function(x) x@effshare)),5) != 1 
             & object[[fl]]@effort > 0) )
      stop("Check fleet '", fl,"', as sum of all effort shares by metier must be 1 if effort > 0." )
    
    
    for (mt in names(object[[fl]]@metiers)) for (st in names(object[[fl]]@metiers[[mt]]@catches)) {
      
      # *.sel: landings.sel + discards.sel == 1
      if ( any(round(landings.sel(object[[fl]]@metiers[[mt]]@catches[[st]]) + 
                     discards.sel(object[[fl]]@metiers[[mt]]@catches[[st]]),5) != 1) ) 
        stop("Check landings.sel and discards.sel for fleet '", fl,"', metier '", mt,"' and stock '", st, "', 
             as their sum must be 1." )
      
      # landings.wt: wage in the landings >= 0
      if ( any(is.na(landings.wt(object[[fl]]@metiers[[mt]]@catches[[st]]))) )
        stop("Required values for landings.wt in fleet '", fl,"', metier '", mt,"' and stock '", st, "'. 
              If not available, set them equal to the mean weights in the stock." )
      
      if ( any(landings.wt(object[[fl]]@metiers[[mt]]@catches[[st]]) < 0) )
        stop("Check landings.wt values for fleet '", fl,"', metier '", mt,"' and stock '", st, "', as they must be >= 0." )
      
      # discards.wt: wage in the discards >= 0
      if ( any(is.na(discards.wt(object[[fl]]@metiers[[mt]]@catches[[st]]))) )
        stop("Required values for discards.wt in fleet '", fl,"', metier '", mt,"' and stock '", st, "'. 
              If not available, set them equal to the mean weights in the landings." )
      
      if ( any(discards.wt(object[[fl]]@metiers[[mt]]@catches[[st]]) < 0) )
        stop("Check discards.wt values for fleet '", fl,"', metier '", mt,"' and stock '", st, "', as they must be >= 0." )
      
    }
    
  }
  
  return(TRUE)
  
}


# data(one)
# checkFleets(oneFl)
# data(oneIt)
# checkFleets(oneItFl)
# data(multi)
# checkFleets(multiFl)
# checkFleets(multiFl, ctrl = multiFlC)
# 
# # Errors
# obj1 <- oneFl
# obj1$fl1@capacity[] <- NA # capacity = NA & no info on fleets.
# checkFleets(obj1)         # outputs a warning
# obj2 <- multiFl
# obj2$fl1@effort
# obj2$fl1@metiers$met1@effshare[,ac(1990),] <- NA # sum != 1, but effort = 0 
# checkFleets(obj2)
# obj2$fl1@metiers$met1@effshare[,ac(1999),] <- 5  # sum != 1, and effort > 0 
# checkFleets(obj2)
# obj3 <- oneFl
# obj3$fl1@metiers$met1@catches$stk1@landings.sel[] <- obj3$fl1@metiers$met1@catches$stk1@discards.sel[] <- 0
# checkFleets(obj3)
# obj4 <- oneFl
# obj4$fl1@metiers$met1@catches$stk1@landings.wt[,5,] <- NA
# checkFleets(obj4)
# obj4$fl1@metiers$met1@catches$stk1@landings.wt[,5,] <- -0.7
# checkFleets(obj4)
# obj5 <- oneFl
# obj5$fl1@metiers$met1@catches$stk1@discards.wt[,5,] <- NA
# checkFleets(obj5)
# obj5$fl1@metiers$met1@catches$stk1@discards.wt[,5,] <- -0.1
# checkFleets(obj5)



#-------------------------------------------------------------------------------
#  checkSRs(object) - checking SRs argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkFLBEIAData
#' @aliases checkSRs

checkSRs <- function(object) {
  
  for (st in names(object)) {
    
    # proportion: 0 < rec distrib along seasons < 1
    if ( any(object[[st]]@proportion < 0) | any(object[[st]]@proportion > 1) ) 
      stop("Check recruitment distribution along seasons for stock '", st,"' (required 0 < proportion < 1)." )
    
    # proportion: sum(recruitment proportion by season)=1
    if ( any(round(seasonSums(object[[st]]@proportion),5) != 1) ) 
      stop("Check recruitment distribution along seasons (i.e. proportion) for stock '", st,"', as their sum must be 1." )
    
    # check: uncertainty > 0
    if ( any(object[[st]]@uncertainty <= 0) ) 
      stop("Check uncertainty in recruitment for stock '", st,"' (required values > 0)." )
    
  }
  
  return(TRUE)
  
}


# data(one)
# checkSRs(oneSR)
# data(oneIt)
# checkSRs(oneItSR)
# data(multi)
# checkSRs(multiSR)
# 
# # Errors
# obj1 <- obj2 <- obj3 <- multiSR
# obj1$stk1@proportion[,,,1,] <- -1000 # proportions > 0
# checkSRs(obj1)
# obj1$stk1@proportion[,,,1,] <- 1000  # proportions < 1
# checkSRs(obj1)
# obj2$stk1@proportion[,,,1:4,] <- 0.5 # sum proportions = 1
# checkSRs(obj2)
# obj3$stk1@uncertainty[1,1,,1,] <- -0.5 # uncertainty> 0
# checkSRs(obj3)



#-------------------------------------------------------------------------------
#  checkBDs(object) - checking BDs argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkFLBEIAData
#' @aliases checkBDs

checkBDs <- function(object) {
  
  for (st in names(object)) {
    
    if(object[[st]]@model=="PellaTom"){
      
      p <- object[[st]]@params["p",,,]
      r <- object[[st]]@params["r",,,]
      K <- object[[st]]@params["K",,,]
      
      if ( any(object[[st]]@alpha<1) || any(as.vector(object[[st]]@alpha) > as.vector(((p/r+1)^(1/p)))) )
        stop("Check alpha parameter for biomass dynamics of stock '", st,"' (required (p/r+1)^(1/p) < alpha < 1)." )
        
    }
  }
  
  return(TRUE)
  
}


# oneBD <- NULL
# checkBDs(oneBD)
# data(multi)
# checkBDs(multiBD)
# 
# # Errors
# obj1 <- obj2 <- obj3 <- multiBD
# obj1$stk2@alpha[1,1,] <- 10 # alpha < 1
# checkBDs(obj1)
# obj2$stk2@alpha[1,1,] <- 
#   (obj2$stk2@params["p",1,1,]/obj2$stk2@params["r",1,1,]+1)^(1/obj2$stk2@params["p",1,1,]) - 1 # alpha > (p/r+1)^(1/p)
# checkBDs(obj2)



#-------------------------------------------------------------------------------
#  checkAdvice(object) - checking advice argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkFLBEIAData
#' @aliases checkAdvice

checkAdvice <- function(object) {
  
  for (st in names(object$quota.share)) {
    
    # quota.shares: sum of all fleets' quota shares == 1
    if ( any(round(quantSums(object$quota.share[[st]]),5) != 1) ) 
      stop("Sum of all quota shares for stock '", st,"' must be 1.")
    
  }
  
  # Quant must be "stock": if not advSum is not working
  if (quant(object$TAC) != "stock")
    stop("quant(object$TAC) must be  'stock', instead of '", quant(object$TAC),"'.")
  
  return(TRUE)
  
}


# data(one)
# checkAdvice(oneAdv)
# data(oneIt)
# checkAdvice(oneItAdv)
# data(multi)
# checkAdvice(multiAdv)
# 
# # Errors
# obj1 <- multiAdv
# obj1$quota.share$stk1[,1,] <- 2 # sum quota shares = 1
# checkAdvice(obj1)



#-------------------------------------------------------------------------------
#  checkObsctrl(object) - checking obs.ctrl argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkFLBEIAData
#' @aliases checkObsctrl

checkObsctrl <- function(object) {
  
  for (st in names(object)) {
    
    sls <- grep(".error",names(object[[st]]$stkObs), value=TRUE)
    
    for (sl in sls)
      if (sl == "ages.error") { # ages.error: sum along ages must be 1
        
        
        if(!is.null(object[[st]]$stkObs$ages.error))
          if ( any(round(apply(object[[st]]$stkObs$ages.error, c(1,3:4), sum),4) != 1) )
            stop("Sum along 2nd dimension of object[['",st,"']]$stkObs$ages.error must be 1.")
        
      } else {                  # *.error: errors > 0
        
        if ( any(object[[st]]$stkObs[[sl]] <= 0) ) 
          stop("Check '", sl,"' in obs.ctrl for stock '", st,"' (required values > 0). 
               For simulating no error, set value to 1." )
        
      }
    
    for (id in names(object[[st]]$indObs)) {
      
      sls <- grep(".error",names(object[[st]]$indObs[[id]]), value=TRUE)
      
      for (sl in sls)
        if (sl == "ages.error") { # ages.error: sum along ages must be 1
          
          
          if(!is.null(object[[st]]$indObs[[id]]$ages.error))
            if ( any(round(apply(object[[st]]$indObs[[id]]$ages.error, c(1,3:4), sum),4) != 1) )
              stop("Sum along 2nd dimension of object[['",st,"']]$indObs[['", id,"']]$ages.error must be 1.")
          
        } else {                  # *.error: errors > 0
          
          if ( any(object[[st]]$indObs[[id]][[sl]] <= 0) ) 
            stop("Check '", sl,"' in obs.ctrl for stock '", st,"' and index '", id,"' (required values > 0). 
                 For simulating no error, set value to 1." )
          
        }
      
    }
    
  }
  
  return(TRUE)
  
}


# data(one)
# checkObsctrl(oneObsC)
# checkObsctrl(oneObsCIndAge)
# checkObsctrl(oneObsCIndBio)
# data(oneIt)
# checkObsctrl(oneItObsC)
# checkObsctrl(oneItObsCIndAge)
# checkObsctrl(oneItObsCIndBio)
# data(multi)
# checkObsctrl(multiObsC)
# 
# # Errors
# obj1 <- oneObsCIndBio
# obj1$stk1$stkObs$land.bio.error[,1,] <- -0.7 # error < 0
# checkObsctrl(obj1)
# obj2 <- oneObsCIndAge
# obj2$stk1$stkObs$ages.error[1,,,] <- 2 # sum ages.error by age = 1
# checkObsctrl(obj2)

