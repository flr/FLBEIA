# check that all flq-s differ only in first (quant) dimension.
# lflq: a list of FLQuant
# dims: the dims we want to check
equal.flq.Dimnames <- function(lflq, dims = 1:6){
    
    dim.flqs <- vector('list', length(lflq))
    names(dim.flqs) <- lflq

    
    for(i in 1:length(lflq)) dim.flqs[[i]] <- dimnames(lflq[[i]])[dims]
    

    for(j in 2:length(lflq)) if(!identical(dim.flqs[[1]], dim.flqs[[j]]))  return(FALSE)
    
    
    return(TRUE)
}
    

# FALSE <- stop("All the input 'FLquant's must share 'year', 'unit', 'season', 'area' and 'iter' dimensions.")

# TEST:
# flq.ank <- FLQuant(dimnames = list(age = 1:2, year = 1999:2002, season = 1:4, iter = 1:5))
# flq.meg <- FLQuant(dimnames = list(age = 0:2, year = 1999:2002, season = 1:4, iter = 1:5))
# flq <- FLQuant(dimnames = list(quant = 'all', year = 1999:2002, season = 1:4, iter = 1:2))
# equal.flq.Dimnames(list(flq.meg, flq.ank), c(2,3)) 
# equal.flq.Dimnames(list(flq.meg, flq.ank), c(1,3)) 
# equal.flq.Dimnames(list(flq.meg, flq.ank, flq), c(1,3)) 
# equal.flq.Dimnames(list(flq.meg, flq.ank, flq), c(2,3)) 
#

