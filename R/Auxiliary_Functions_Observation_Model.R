#' @rdname landStock.f
#' @aliases tlandStock
#Landings

#tlandStock
#Total landings of a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock

tlandStock <- function(obj, stknm){
    aux <- 0
    for(f in obj){
        for(m in f@metiers){
            if(!(stknm %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stknm]]@landings
                res[is.na(res)] <- 0
                next
            }
            resf <- m@catches[[stknm]]@landings
            resf[is.na(resf)] <- 0
            res <- res + resf
        }
    }
    return(res)
}

#landStock()
#landings in numbers at age across fleets and metiers
#function landStock already exists

#wtalStock
#Mean weight at age in landings for a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock

#' Mean weight at age in landings for a stock across fleets and metiers
#' 
#' @param obj An object of class FlFleetsExt.
#' @param stknm Character. The name of the stock for which we want to calculate mean weight-at-age.
#
#' @return A FLQuant with mean weight-at-age values. 
#' 
#' @seealso \code{\link{wtadStock}} 
#' 

wtalStock <- function(obj, stknm)
    {
    aux <- 0
    cnt <- 1
    for(f in obj)
        {
        for(m in f@metiers)
           {
           if(!(stknm %in% catchNames(m)))
             next
           if(aux == 0)
             {
             aux <- 1
             res <- m@catches[[stknm]]@landings.wt
             res[is.na(res)] <- 0
             next
             }
           cnt <- cnt + 1
           resf <- m@catches[[stknm]]@landings.wt
           resf[is.na(resf)] <- 0
           res <- res + resf
           }
        }
    res <- res/cnt # This is a normal mean ant it should be weighted!!
    return(res)
}

#tdiscStock
#Total discards of a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock
#' @rdname landStock.f
#' @aliases tdiscStock
tdiscStock <- function(obj, stknm){
    aux <- 0
    for(f in obj){
        for(m in f@metiers){
            if(!(stknm %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stknm]]@discards
                res[is.na(res)] <- 0
                next
            }
            resf <- m@catches[[stknm]]@discards 
            resf[is.na(resf)] <- 0
            res <- res + resf
        }
    }
    return(res)
}

#discnStock
#Discards in numbers at age of a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock
discnStock <- function(obj, stknm){
    aux <- 0
    for(f in obj){
        for(m in f@metiers){
            if(!(stknm %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stknm]]@discards.n
                res[is.na(res)] <- 0
                next
            }
            resf <- m@catches[[stknm]]@discards.n
            resf[is.na(resf)] <- 0
            res <- res + resf
        }
    }
    return(res)
}

#wtadStock
#Mean weight at age in landings for a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock

#' Mean weight at age in discards for a stock across fleets and metiers
#' 
#' @param obj An object of class FlFleetsExt.
#' @param stknm Character. The name of the stock for which we want to calculate mean weight-at-age.
#
#' @return A FLQuant with mean weight-at-age values.
#' 
#' @seealso \code{\link{wtalStock}} 
#' 

wtadStock <- function(obj, stknm)
    {
    aux <- 0
    cnt <- 1
    for(f in obj)
        {
        for(m in f@metiers)
           {
           if(!(stknm %in% catchNames(m)))
             next
           if(aux == 0)
             {
             aux <- 1
             res <- m@catches[[stknm]]@discards.wt
             res[is.na(res)] <- 0
             next
             }
           cnt <- cnt + 1
           resf <- m@catches[[stknm]]@discards.wt 
           resf[is.na(resf)] <- 0
           res <- res + resf
           }
        }
    res <- res/cnt
    return(res)
}
