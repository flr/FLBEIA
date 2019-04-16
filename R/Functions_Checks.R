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



#' @title Checking arguments of FLBEIA function
#' @description This functions checks wether some conditions are met for the inputs 
#'              in the different arguments of the FLBEIA function.
#'              Specifically there are checking functions for: 
#'              biols, fleets, SRs, advice, and obs.ctrl arguments.
#' 
#' Functions available:
#' \itemize{
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
#' @param object An object of the appropiate class depending on the FLBEIA argument to be checked:
#'               - biols   : A FLBiols object.
#'               - fleets  : A FLFleetsExt object.
#'               - SRs     : A list of FLSRsim objects.
#'               - advice  : A list with two FLQuant elements, TAC and quota.share.
#'               - obs.ctrl: A list.
#'
#' @return An error message if there is any error detected.
#' 


#-------------------------------------------------------------------------------
#  checkBiols(object) - checking biols argument of FLBEIA function
#-------------------------------------------------------------------------------

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
#' @rdname checkBiols
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
# obj1$fl1@capacity[] <- NA # capacity = NA
# checkFleets(obj1)
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
#' @rdname checkBiols
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
#  checkAdvice(object) - checking advice argument of FLBEIA function
#-------------------------------------------------------------------------------
#' @rdname checkBiols
#' @aliases checkAdvice

checkAdvice <- function(object) {
  
  for (st in names(object$quota.share)) {
    
    # quota.shares: sum of all fleets' quota shares == 1
    if ( any(round(quantSums(object$quota.share[[st]]),5) != 1) ) 
      stop("Sum of all quota shares for stock '", st,"' must be 1.")
    
  }
  
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
#' @rdname checkBiols
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

