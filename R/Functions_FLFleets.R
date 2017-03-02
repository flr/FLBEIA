#-------------------------------------------------------------------------------
#              - Functions to ease the work with FLFLeets -
#
#    - landStock: Extract total landings at age of a stock from a FLFleets obj.
#    - discStock: Extract total discards at age of a stock from a FLFleets obj.
#    - catchStock: Extract total catch at age of a stock from a FLFleets obj.
#    - landWStock, discWStock, catchWStock: The same but in weight.
#    - stock.fleetInfo: Return a matrix [ns,fl_mt]with rownames equal to the 
#           name of the stocks and colnames equal to fleet names '&&' metier names.
#           if element (i,j) = 0, ith stock is not caught by fleet/metier flname&&metiername.
#           
#
# Created: 19/10/2010 16:15:49
# Author: Dorleta Garcia
# Changed: 12/01/2011 15:45:27
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#  landStock(obj, stock)
#  obj : FLFleets object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
landStock <- function(obj, stock){
    
    aux <- 0
    res <- FLQuant()
    
    for(f in obj){
        for(m in f@metiers){

            if(!(stock %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stock]]@landings.n 
                res[is.na(res)] <- 0
                next
            }
            
            resf <- m@catches[[stock]]@landings.n 
            resf[is.na(resf)] <- 0
            res <- res + resf
        }
    }
    return(res)
}
        
#-------------------------------------------------------------------------------
#  discStock(obj, stock)
#  obj : FLFleets object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
discStock <- function(obj, stock){
    aux <- 0
    res <- FLQuant()
    
    for(f in obj){
        for(m in f@metiers){

            if(!(stock %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stock]]@discards.n 
                res[is.na(res)] <- 0
                next
            }   
            resf <- m@catches[[stock]]@discards.n 
            resf[is.na(resf)] <- 0
            res <- res + resf
            
        }
    }
    return(res)

}

#-------------------------------------------------------------------------------
#  catchStock(obj, stock)
#  obj : FLFleets object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
catchStock <- function(obj, stock){
    return(discStock(obj, stock) + landStock(obj,stock))
}


#-------------------------------------------------------------------------------
#  landWStock(obj, stock)
#  obj : FLFleets object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
landWStock <- function(obj, stock){
    
    aux <- 0
    res <- FLQuant()
    
    for(f in obj){
        for(m in f@metiers){

            if(!(stock %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stock]]@landings.n*m@catches[[stock]]@landings.wt
                res[is.na(res)] <- 0
                next
            }
            resf <- m@catches[[stock]]@landings.n*m@catches[[stock]]@landings.wt 
            resf[is.na(resf)] <- 0
            res <- res + resf
            
        }
    }
    return(res)
}
        
#-------------------------------------------------------------------------------
#  discWStock(obj, stock)
#  obj : FLFleets object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
discWStock <- function(obj, stock){
    aux <- 0
    res <- FLQuant()
    
    for(f in obj){
        for(m in f@metiers){

            if(!(stock %in% catchNames(m)))
                next
            if(aux == 0){
                aux <- 1
                res <- m@catches[[stock]]@discards.n*m@catches[[stock]]@discards.wt 
                res[is.na(res)] <- 0
                next
            }
            resf <- m@catches[[stock]]@discards.n*m@catches[[stock]]@discards.wt 
            resf[is.na(resf)] <- 0
            res <- res + resf
            
        }
    }
    return(res)

}

#-------------------------------------------------------------------------------
#  catchWStock(obj, stock)
#  obj : FLFleets object. 
#  Output: An FLQuant(age, year, unit, n_season,1,iter) with the total catch of 
#      the 'stock' in the FLFleets objetc.
#-------------------------------------------------------------------------------
#' @rdname landStock.f
catchWStock <- function(obj, stock){
    return(discWStock(obj, stock) + landWStock(obj,stock))
}


#-------------------------------------------------------------------------------
#    - stock.fleetInfo: Return a matrix [ns,fl_mt]with rownames equal to the 
#           name of the stocks and colnames equal to fleet names '&&' metier names.
#           if element (i,j) = 0, ith stock is not caught by fleet/metier flname&&metiername.
#-------------------------------------------------------------------------------
#' stock.fleetInfo
#'
#' Indicates which stocks are caugth by each fleet-metier combination.
#'
#' @return Return a matrix with rownames equal to the stocks names and colnames equal to names of fleet and metier.
#' If element (i,j) is equal to 0, the stock (i) is not caught by fleet/metier (j). 
#'
#' @param fleets is an objetc of class FLFleetsExt.
#'
#' @examples
#' \dontrun{
#' data(multi)
#' stock.fleetInfo(fl1)
#' }



stock.fleetInfo <- function(fleets){
   # sts  <- catchNames(fleets)
    
    stfl <- lapply(fleets, catchNames)
    sts <- NULL 
    for(ss in stfl) sts <- c(sts, ss)
    sts <- unique(sts)

    nf <- length(fleets)
    flmt <- unlist(lapply(1:nf, function(x) paste(fleets[[x]]@name, names(fleets[[x]]@metiers), sep = '&&')))
    res <- matrix(0, length(sts), length(flmt), dimnames = list(stock = sts, fleet_met = flmt))

    for(fm in flmt){
   # print(fm)
        xx   <- strsplit(fm, "&&")[[1]]
        fmst <- names(fleets[[xx[1]]]@metiers[[xx[2]]]@catches)
        res[fmst,fm] <- 1
    }
    return(res)
}





