#-------------------------------------------------------------------------------
#   propagateFLF: Propagate an FLFleet object
#       object: An FLFleet object.
#
#   fillIterFLF: Fill iter-s in FLFleet oobject
#       
# Dorleta Garcia
# created: 03/11/2011 08:34:38
#-------------------------------------------------------------------------------

propagateFLF <- function(object, iter, fill.iter){

    nmet <- length(object@metiers)
    
    catches <- lapply(object@metiers, function(x) 
                        lapply(x@catches, propagate, iter = iter, fill.iter = fill.iter))
   
    metiers <- vector('list', nmet)
     
    for(mt in 1:nmet){
                   
        metiers[[mt]] <- FLMetierExt(name = object[[mt]]@name, desc = object[[mt]]@desc, range = object[[mt]]@range, gear = object[[mt]]@gear,
                                 effshare = propagate(object@metiers[[mt]]@effshare, iter = iter, fill.iter = fill.iter),
                                    vcost = propagate(object@metiers[[mt]]@vcost, iter = iter, fill.iter = fill.iter),
                                  catches = FLCatchesExt(catches[[mt]]))  
    }
    
    names(metiers) <- names(object@metiers)
                                
    fleet <- FLFleetExt(name = object@name, desc = object@desc, range = object@range,
                      effort = propagate(object@effort, iter = iter, fill.iter = fill.iter), 
                       fcost = propagate(object@fcost, iter = iter, fill.iter = fill.iter),
                    capacity = propagate(object@capacity, iter = iter, fill.iter = fill.iter),
                   crewshare = propagate(object@crewshare, iter = iter, fill.iter = fill.iter),
                     metiers = FLMetiersExt(metiers))
    
    return(fleet)
}
                      
       


# fillIter

fillIterFLF <- function(obj1, obj2, iter){

    mnms <- names(obj1@metiers)
    
    for(mt in mnms){
        
        cnms <- names(obj1[[mt]]@catches)
         
        for(ct in cnms){
            iter(obj1@metiers[[mt]]@catches[[ct]],iter) <- obj2@metiers[[mt]]@catches[[ct]] 
        }
        
        iter(obj1@metiers[[mt]]@vcost,iter)    <- obj2@metiers[[mt]]@vcost
        iter(obj1@metiers[[mt]]@effshare,iter) <- obj2@metiers[[mt]]@effshare 
    }

    iter(obj1@effort,iter)    <- obj2@effort 
    iter(obj1@fcost,iter)     <- obj2@fcost 
    iter(obj1@capacity,iter)  <- obj2@capacity 
    iter(obj1@crewshare,iter) <- obj2@crewshare   
    
    return(obj1)
}




