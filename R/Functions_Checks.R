#-------------------------------------------------------------------------------
#              - Functions to run basic checks.    
#
# Created: 17/01/2011 14:47:19
# Author: Dorleta Garcia
# Changed: 17/01/2011 14:47:24
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

