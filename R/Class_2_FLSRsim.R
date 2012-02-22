# An object to simulate recruitment.

validFLSRsim <- function(object){

   ## All FLQuant objects must have same dimensions, including 'iter'.
   Dim <- dim(object@rec)

   s.  <-list("ssb","uncertainty", "proportion")

   for (i. in s.)	{
	   t. <- slot(object,i.)
	   if (is.FLQuant(t.) & !all(dim(t.) == Dim))
	      stop(paste("FLQuant dimensions wrong for ", i.))
	   }
  
  ## FLQuant in covar all the same dimensions.
  for(i in names(object@covar)){
	   if (!all(dim(object@covar[[i]]) == Dim))
	      stop(paste("FLQuant dimensions in 'covar' list wrong for ", i.))   
  }

  # Params.
   if(!(all(dim(object@params)[-1] == Dim[c(2,4,6)])))  
          stop(cat("Wrong dimension in 'params', it should be equalt to : npar", Dim[c(2,4,6)], '\n'))   
  
  # timelag
    if(!(all(dim(object@timelag) == c(2,Dim[4]))))  
          stop(cat("Wrong dimension in 'timelag', it should be equal to:",c(2,Dim[4]), '\n'))  
    
  # model
    if(length(object@model) != 1)  
          stop("Wrong length in 'model', it must be equal  1")  
    if(length(grep('~', object@model)) == 0 &  class(eval(call(object@model))[[2]]) != 'formula')
        stop("Specified 'model' is not defined in 'FLCore' or the specified formula is not correct")

  # season in timelag must be >= 1 and <= ns
  if(any(object@timelag[2,] < 1 | object@timelag[2,] > Dim[4]))
        stop(cat('season in timelag must be between 1 and', Dim[4], '\n')) 

   # Everything is fine
   return(TRUE)
   }
   
   
setClass("FLSRsim",
	representation(
		"FLComp",
        rec               = "FLQuant",        # [1,ny,1,ns,1,it]
		ssb               = "FLQuant",        # [1,ny,1,ns,1,it]
		covar             = "FLQuants",       # [1?¿,ny,1,ns,1,it]
		uncertainty       = "FLQuant",        # [1,ny,1,ns,1,it]
		proportion        = "FLQuant",        # [1,ny,1,ns,1,it]
		model             = "character", # [it] - different model by iteration.
		params            = "array",          # array[param, year, season, iteration]    # year in order to model regime shifts.
		timelag           = "matrix"         # [2,ns]
      ),
	prototype=prototype(
		name     =character(0),
		desc     =character(0),
		range    =unlist(list(min=NA, max=NA, plusgroup=NA, minyear=1, maxyear=1)),
        rec               = FLQuant(),        # [1,ny,1,ns,1,it]
		ssb               = FLQuant(),        # [1,ny,1,ns,1,it]
		covar             = FLQuants(),       # [1?¿,ny,1,ns,1,it]
		uncertainty       = FLQuant(),        # [1,ny,1,ns,1,it]
		proportion        = FLQuant(),        # [1,ny,1,ns,1,it]
		model             = as.character(NA), # [it] - different model by iteration.
		params            = array(),          # array[param, year, season, iteration]    # year in order to model regime shifts.
		timelag           = matrix(c(0,1),2,1, dimnames = list(c('year', 'season'),season = 1)),         # [2,ns]  year = how many years before, season = which season.
	    validity=validFLSRsim
))

setValidity("FLSRsim", validFLSRsim)
remove(validFLSRsim)	# We do not need this function any more
#invisible(createFLAccesors("FLSRSim", exclude=c('name', 'desc', 'range'))) # }}}

FLSRsim <- function(...){
    a <- new('FLSRsim')
    x <- list(...)
    attach(x, pos = 2, warn.conflicts = FALSE)
    slots <- names(x)
    
    if(all(slots %in% c('name', 'desc', 'range'))){   # => NO DIMENSIONS
        for(sl in slots){
            slot(a,sl) <- x[sl][[1]]
        }
        return(a)
    }
    
    # ELSE, CORRECT DIMENSIONS ARE NEEDED IN THE DIMENSIONAL SLOTS!!!
    ny <- ns <- ni <- 1
    nmy  <- nmi <- 1
    nms <- 'all'
    nmparams <- 'a'
     
    quants <- c('rec', 'ssb', 'uncertainty', 'proportion')
    
    if(any(slots %in% quants)){  
        Dim    <- dim(get(slots[which(slots  %in% quants)[1]], pos = 2))
        Dimnms <- dimnames(get(slots[which(slots  %in%quants)[1]], pos = 2))
        ny <- Dim[2]; nmy <- Dimnms[[2]]
        ns <- Dim[4]; nms <- Dimnms[[4]]
        ni <- Dim[6]; nmi <- Dimnms[[6]]
    }
    else{  # slots in  "covar"    "model"       "params"      "timelag"     "name"        "desc"   "range" 
        if(any(slots == 'covar')){
            Dim <- dim(x[['covar']][[1]])
            Dimnms <- dimnames(x[['covar']][[1]])
            ny <- Dim[2]; nmy <- Dimnms[[2]]
            ns <- Dim[4]; nms <- Dimnms[[4]]
            ni <- Dim[6]; nmi <- Dimnms[[6]]
        }
        else{
            if(any(slots == 'params')){
                Dim    <- dim(x[['params']])
                Dimnms <- dimnames(x[['params']])
                ny <- Dim[2]; nmy <- Dimnms[[2]]
                ns <- Dim[3]; nms <- Dimnms[[3]]
                ni <- Dim[4]; nmi <- Dimnms[[4]]
            }
            else{
                if(any(slots == 'timelag')){
                    ns <- dim(x[['timelag']])[2]; nms <- dimnames(x[['timelag']])[[2]]
                }
                else{ # => model
                    Dim <- length(x[['model']])
                    ni <- Dim ; nmi <- 1:Dim
                }            
            }
        }
    }      
        
    if('model' %in% slots){
        nmparams <- all.vars(get(x[['model']], pos = 2)()[[2]])
        nmparams <- nmparams[!(nmparams %in% c('rec', 'ssb'))]
    }
    
    detach()
    
    for(sl in slots){ # fill the slots that are in the call.
        slot(a,sl) <- x[sl][[1]]
    }
    
    # dimensional slots that are not in the call.  
    # if covar is not in the call => covar :: EMPTY.
    noin <- slotNames(a)[-which((slotNames(a) %in% slots))]
    noin <- noin[! noin %in% c('name', 'desc', 'range', 'covar')]
   
    for(sl in noin){
        if(sl %in% quants){
             y <- ifelse(sl %in% c('proportion', 'uncertainty'), 1, NA) 
            xx <- ifelse(sl == 'rec', 0, 'all')
            slot(a,sl) <- FLQuant(y,dim = c(1,ny,1,ns,1,ni), dimnames = list(age = 'all', year = nmy, season = nms, iter = nmi))  
        }   
        else
            if(sl == 'model') slot(a,sl) <- 'geomean'
                else
                    if(sl == 'params') 
                        slot(a, sl) <- array(dim = c(length(nmparams),ny,ns,ni), dimnames = list(params = nmparams, year = nmy, season = nms, iter = nmi))
                    else # => timelag
                        slot(a,sl) <- matrix(c(rep(0,ns),1:ns),2,ns, dimnames = list(c('year', 'season'), nms), byrow = TRUE)
   }

    validObject(a)
    return(a)
}


# SRsim   {{{
# It fills in  the recruitment slot for 
# Only 1 year and 1 season. Otherwise we would have ny*ns combinations and for simulation 
# we simulate one by one! 
# iter >= 1.

SRsim <- function(object, year = 1, season = 1, iter = 'all')  # year and season either numeric (position) or character (name)
  {

    Dim    <- dim(object@rec)
    dimnms <- dimnames(object@rec)
    
    # If year/season/iter numerics => indicate position 
    # else names => get positions.
    
    if(length(year) > 1 | length(season) > 1)
    stop('Only one year and season is allowed' )
    
    # 'year' dimension.
    yr <- year
    if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
    if(length(yr) == 0) stop('The year is outside object time range')  
    
    # 'season' dimension.
    ss <- season
    if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
    if(length(ss) == 0) stop('The season is outside object season range')  

     # 'iter' dimension.
    if(iter == 'all') it <- 1:Dim[6]
    else if(is.character(iter)) it <- which(dimnms[[6]] %in% iter)
    if(length(it) == 0 | length(it) < length(iter)) stop('Some of all iterations are outside object iteration range')  

    # model call
    if(length(grep('~', object@model)) == 0)
        model <- eval(call(object@model))[[2]]
    else # character but 'formula' 
        model <- formula(object@model)
    
    model <- as.list(model)[[3]]
     
    # Extract SSB season and year according to timelag, 
    # Extract two recruitments, recruitment in previous season and recruitment 
    # exactly one year before.
    ns  <- dim(object@ssb)[4]
    ss1 <- ifelse(ns == 1, 1, ifelse(ss == 1, ns, ss-1))
    yr1 <- ifelse(ss > 1, yr, yr-1)
    
    datam <- list(ssb       = c(object@ssb[,yr - object@timelag[1,ss],, object@timelag[2,ss],]),
                  rec.prevS = c(object@rec[,yr1,,ss1]),
                  rec.prevY = c(object@rec[,yr-1,,ss]))

    # Extract covars.
    for(i in names(object@covar))
        datam[[i]] <-  c(object@covar[[i]][,yr,, ss,])
    
    # Extract params
    for(i in dimnames(object@params)[[1]])
        datam[[i]] <-  c(object@params[i,yr,ss,])

    res <- c(eval(model, datam))   # valid for 1 year, 1 season and 'N' iterations
    
    object@rec[,yr,,ss,] <- res*object@proportion[,yr,,ss,]*object@uncertainty[,yr,,ss,]
    
    return(object)
    
}   # }}}


## iter {{{
setMethod("iter", signature(object="FLSRsim"),
	  function(object, iter)
	  {
		# FLQuant slots
		names <- names(getSlots(class(object))[getSlots(class(object))=="FLQuant"])
		for(s in names) 
		{
			if(dims(slot(object, s))$iter == 1)
				slot(object, s) <- iter(slot(object, s), 1)
			else
				slot(object, s) <- iter(slot(object, s), iter)
		}
		# covar
		if(length(object@covar) > 0) slot(object, 'covar') <- iter(slot(object, 'covar'), iter)
        
        #params
        slot(object, 'params') <-  slot(object, 'params')[,,,iter,drop=F]
        dimnames(slot(object, 'params'))[[4]] <- 1:length(iter) 
         
		return(object)
	  }
) # }}}


