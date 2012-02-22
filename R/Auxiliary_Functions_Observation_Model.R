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
            res <- res + m@catches[[stknm]]@landings
            res[is.na(res)] <- 0
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
           res <- res + m@catches[[stknm]]@landings.wt
           res[is.na(res)] <- 0
           }
        }
    res <- res/cnt
    return(res)
}

#tdiscStock
#Total discards of a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock
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
            res <- res + m@catches[[stknm]]@discards
            res[is.na(res)] <- 0
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
            res <- res + m@catches[[stknm]]@discards.n
            res[is.na(res)] <- 0
        }
    }
    return(res)
}

#wtadStock
#Mean weight at age in landings for a stock across fleets and metiers
#obj : an object of class FlFleetsExt
#stknm : character, the name of the stock
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
           res <- res + m@catches[[stknm]]@discards.wt
           res[is.na(res)] <- 0
           }
        }
    res <- res/cnt
    return(res)
}
