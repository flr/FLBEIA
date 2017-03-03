#-------------------------------------------------------------------------------
#
#' Extract landings, discard or catch at age 
#' 
#' Extract landings, discards or catch at age of a stock from a FLFleetExt or FLFleetsExt object.
#'
#' @details
#'\itemize{
#'      \item{landStock.f:} Extract total landings at age of a stock from a FLFleetExt obj.
#'      \item{discStock.f:} Extract total discards at age of a stock from a FLFleetExt obj. 
#'      \item{catchStock.f:} Extract total catch at age of a stock from a FLFleetExt obj.    
#'      \item{landWStock.f:} Extract total landings at age in weight of a stock from a FLFleetExt obj.
#'      \item{discWStock.f:} Extract total discards at age in weight of a stock from a FLFleetExt obj.
#'      \item{catchWStock.f:} Extract total catch at age in weight of a stock from a FLFleetExt obj. 
#'      \item{landStock:} Extract total landings at age of a stock from a FLFleetExts obj.
#'      \item{discStock:} Extract total discards at age of a stock from a FLFleetExts obj. 
#'      \item{catchStock:} Extract total catch at age of a stock from a FLFleetExts obj.    
#'      \item{landWStock:} Extract total landings at age in weight of a stock from a FLFleetExts obj.
#'      \item{discWStock:} Extract total discards at age in weight of a stock from a FLFleetExts obj.
#'      \item{catchWStock:} Extract total catch at age in weight of a stock from a FLFleetExts obj.                       
#'      \item{tlandStock:} Total landings of a stock across fleets and metiers                
#'      \item{tdiscStock:} Total discard of a stock across fleets and metiers  
#'}


#' @param obj a FLFleetExt (.f) or FLFleetsExt object 
#' @param stock stock name
#' @param stknm stock name 
#
#' @return A FLQuant object with landings, discards or catch (total or at age). 

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' data(one)
#' landWStock.f(oneFl$fl1,"stk1")
#' landWStock(oneFl,"stk1")
#' }
#
#              - Functions to ease the work with FLFLeetExt -
#
#    - landStock.f: Extract total landings at age of a stock from a FLFleet obj.
#    - discStock.f: Extract total discards at age of a stock from a FLFleet obj.
#    - catchStock.f: Extract total catch at age of a stock from a FLFleet obj.
#    - landWStock.f, discWStock.f, catchWStock.f: The same but in weight.
#           
#
# Created: 30/11/2010 13:23:35
# Author: Dorleta Garcia
# Changed: 12/01/2011 15:48:43
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#  landStock.f(obj, stock)
#  obj : FLFleet object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------

# @export
landStock.f <- function(obj, stock){
    
    aux <- 0
    res <- FLQuant()
    
    for(m in obj@metiers){

        if(!(stock %in% catchNames(m)))
            next
        if(aux == 0){
            aux <- 1
            res <- m@catches[[stock]]@landings.n 
            res[is.na(res)] <- 0
            next
        }
        res <- res + m@catches[[stock]]@landings.n 
        res[is.na(res)] <- 0
    }
    return(res)
}
        
#-------------------------------------------------------------------------------
#  discStock.f(obj, stock)
#  obj : FLFleet object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------

#' @rdname landStock.f
discStock.f <- function(obj, stock){
    aux <- 0
    res <- FLQuant()
    
    for(m in obj@metiers){
    
        if(!(stock %in% catchNames(m))) next
        
        if(aux == 0){
            aux <- 1
            res <- m@catches[[stock]]@discards.n 
            res[is.na(res)] <- 0
            next
        }
        res <- res + m@catches[[stock]]@discards.n 
        res[is.na(res)] <- 0
    }
    return(res)
}

#-------------------------------------------------------------------------------
#  catchStock.f(obj, stock)
#  obj : FLFleet object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------

#' @rdname landStock.f
catchStock.f <- function(obj, stock){
    return(discStock.f(obj, stock) + landStock.f(obj,stock))
}


#-------------------------------------------------------------------------------
#  landWStock.f(obj, stock)
#  obj : FLFleet object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
landWStock.f <- function(obj, stock){
    
    aux <- 0
    res <- FLQuant()
    
    for(m in obj@metiers){

        if(!(stock %in% catchNames(m))) next
        
        if(aux == 0){
            aux <- 1
            res <- m@catches[[stock]]@landings.n*m@catches[[stock]]@landings.wt
            res[is.na(res)] <- 0
            next
        }
        res <- res + m@catches[[stock]]@landings.n*m@catches[[stock]]@landings.wt 
        res[is.na(res)] <- 0
    }

    return(res)
}
        
#-------------------------------------------------------------------------------
#  discWStock.f(obj, stock)
#  obj : FLFleet object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
discWStock.f <- function(obj, stock){
    aux <- 0
    res <- FLQuant()
    
    for(m in obj@metiers){

        if(!(stock %in% catchNames(m))) next
        
        if(aux == 0){
            aux <- 1
            res <- m@catches[[stock]]@discards.n*m@catches[[stock]]@discards.wt 
            res[is.na(res)] <- 0
            next
        }
        res <- res + m@catches[[stock]]@discards.n*m@catches[[stock]]@discards.wt 
        res[is.na(res)] <- 0
    }
    return(res)
}

#-------------------------------------------------------------------------------
#  catchWStock.f(obj, stock)
#  obj : FLFleet object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
catchWStock.f <- function(obj, stock){
    return(discWStock.f(obj, stock) + landWStock.f(obj,stock))
}






