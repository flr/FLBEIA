#------------------------------------------------------------------------------#
#                           Join Iterations 
#  
#   * joinIter(object, files, directory, Niters)
#
#  
#   object: character. The name of the object that must be 'joined'. The object 
#           _must_ be the output of BEIA function.
#   files: a character vector with the names of the files.
#   directory: the directory where the files are stored, (the default is the 
#             current directory).       
#   Niters: a numeric vector with the number of iterations per object. If length = 1
#          it is assumed that all objects have the same number of iterations.
#   elements: the elements of the objects that must be joined. The default is to 
#           join all the objects.
#   note: The files must contain a single object
# 
#   
# Dorleta García
# Created: 03/11/2011 07:42:23
# Changed: 03/11/2011 07:42:28
#------------------------------------------------------------------------------#

joinIter <- function(object, files, directory = NULL, Niters = 1, elements = 'all', advice.ext = 'TAC', fleets.ctrl.ext = 'seasonal.share'){

    nit <- ifelse(length(Niters) == 1, length(files)*Niters, prod(Niters))

    if(length(Niters) == 1)  Niters <- rep(Niters, length(files))
    
    # load the one object and propagate it to 'nit'.
    if(!is.null(directory)) setwd(directory)
    
    objnam <- object
    
    load(files[1])
    
    object <- get(objnam)
    
    # Propagate the object.
    
    advice <- fleets.ctrl <- list()
    
    biols   <- lapply(object$biols, propagate, nit, fill.iter = FALSE)
    fleets  <- lapply(object$fleets, propagateFLF, nit, fill.iter = FALSE)
    covars  <- lapply(object$covars, propagate, nit, fill.iter = FALSE)
    stocks  <- lapply(object$stocks, propagate, nit, fill.iter = FALSE)
    indices <- lapply(object$indices, function(x) FLIndices(lapply(x, propagate, nit, fill.iter = FALSE)))
    advice[[advice.ext]]            <- propagate(object$advice[[advice.ext]], nit, fill.iter = FALSE)
    fleets.ctrl[[fleets.ctrl.ext]]  <- lapply(object$fleets.ctrl[[fleets.ctrl.ext]], propagate, nit, fill.iter = FALSE)
    
    k <- 2
    iter0 <- Niters[1] + 1

    for(fl in files[-1]){
        iter.sel <- iter0:(iter0 + Niters[k] - 1) 
        load(fl)
        
        print(iter.sel)

        object <- get(objnam)
        
        # biols
        for(i in names(biols)) iter(biols[[i]], iter.sel) <- object$biols[[i]] 
        
        # fleets 
        for(i in names(fleets)) fleets[[i]] <- fillIterFLF(fleets[[i]], object$fleets[[i]],  iter.sel) 
        
        # covars
        for(i in names(covars)) iter(covars[[i]], iter.sel) <- object$covars[[i]]
        
        # stocks 
        for(i in names(stocks)) iter(stocks[[i]], iter.sel) <- object$stocks[[i]]
        
        # indices
        for(i in names(indices)){
            for(j in names(indices[[i]])){ 
                iter(indices[[i]][[j]], iter.sel) <- object$indices[[i]][[j]]
            }
        }
        
        # advice[[advice.ext]]
        iter(advice[[advice.ext]], iter.sel) <- object$advice[[advice.ext]]
        
        # advice[[advice.ext]]
         for(i in names(fleets.ctrl[[fleets.ctrl.ext]]))  
            iter(fleets.ctrl[[fleets.ctrl.ext]][[i]], iter.sel) <- object$fleets.ctrl[[fleets.ctrl.ext]][[i]]
        
        iter0 <- ifelse(length(iter.sel) == 1, iter.sel, iter.sel[2]) + 1
        k <- k+1
    }
    
    fleets <- FLFleetsExt(fleets)
    
    res <- list(biols = biols, fleets = fleets, covars = covars, stocks = stocks, indices = indices, fleets.ctrl = fleets.ctrl, advice = advice)
        
    return(res)    
}
    
     
    
    
    
    
    
    
    

