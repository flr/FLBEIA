#-------------------------------------------------------------------------------  
#  FLBDsim class.
# Created: Sonia Sánchez - 16/08/2010 12:37:18
# Changed: 26/10/2010 08:40:46 (dorleta garcia)
#------------------------------------------------------------------------------- 

# An object to simulate biomass (when using production models).

validFLBDsim <- function(object){

   ## All FLQuant objects must have same dimensions, including 'iter'.
   Dim <- dim(object@biomass)

   s.  <-list("biomass","catch","uncertainty")

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
  
  # model
    if(length(object@model) != 1)  
          stop("Wrong length in 'model', it must be equal  1")  
    if(length(grep('~', object@model)) == 0 &  class(eval(call(object@model))[[2]]) != 'formula')
        stop("Specified 'model' is not defined in 'FLCore' or the specified formula is not correct")

   # check that the number of params and the names are consistent with the selected model
#    npar <- all.vars(get(object@model, pos = 1)()[[2]])
#    npar <- npar[!(npar %in% c('biomass','bprod'))]
#    if ( dim(object@params)[1] != length(npar) )
#      stop(cat("For default model 'PellaTom', the number of parameters should be 3 (", npar,")", '\n'))
#    for (p in npar) if ( !(p %in% dimnames(object@params)$params))
#      stop(cat("Parameters should be:", npar, '\n'))
   
   # Everything is fine
   return(TRUE)

   }
   
   
setClass("FLBDsim",
	representation(
		"FLComp",
    biomass           = "FLQuant",        # [1,ny,1,ns,1,it]
    catch             = "FLQuant",        # [1,ny,1,ns,1,it]
		covar             = "FLQuants",       # [1,ny,1,ns,1,it]
		uncertainty       = "FLQuant",        # [1,ny,1,ns,1,it]
		model             = "character",      # [it] - different model by iteration.
		params            = "array"           # array[param, year, season, iteration]    # year in order to model regime shifts.
  ),
	prototype=prototype(
		name     =character(0),
		desc     =character(0),
		range    =unlist(list(min=NA, max=NA, plusgroup=NA, minyear=1, maxyear=1)),
    biomass           = FLQuant(),        # [1,ny,1,ns,1,it]
		catch             = FLQuant(),        # [1,ny,1,ns,1,it]
		covar             = FLQuants(),       # [1,ny,1,ns,1,it]
		uncertainty       = FLQuant(),        # [1,ny,1,ns,1,it]
		model             = as.character(NA), # [it] - different model by iteration.
		params            = array(),          # array[param, year, season, iteration]    # year in order to model regime shifts.
	    validity=validFLBDsim
))

setValidity("FLBDsim", validFLBDsim)
remove(validFLBDsim)	# We do not need this function any more
#invisible(createFLAccesors("FLBDSim", exclude=c('name', 'desc', 'range'))) # }}}

FLBDsim <- function(...){
    a <- new('FLBDsim')
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
     
    quants <- c('biomass', 'catch', 'uncertainty')
    
    if(any(slots %in% quants)){  
        Dim    <- dim(get(slots[which(slots  %in% quants)[1]], pos = 2))
        Dimnms <- dimnames(get(slots[which(slots  %in%quants)[1]], pos = 2))
        ny <- Dim[2]; nmy <- Dimnms[[2]]
        ns <- Dim[4]; nms <- Dimnms[[4]]
        ni <- Dim[6]; nmi <- Dimnms[[6]]
    }
    else{  # slots in   "covar"   "model"   "params"   "name"   "desc"   "range" 
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
            else{ # => model
                Dim <- length(x[['model']])
                ni <- Dim ; nmi <- 1:Dim
            }            
        }
    }      
        
    if('model' %in% slots){
        nmparams <- all.vars(get(x[['model']], pos = 1)()[[2]])
        nmparams <- nmparams[!(nmparams %in% c('biomass','bprod'))]
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
             y <- ifelse(sl %in% c('uncertainty'), 1, NA) 
            xx <- ifelse(sl == 'biomass', 0, 'all')
            slot(a,sl) <- FLQuant(y,dim = c(1,ny,1,ns,1,ni), dimnames = list(age = 'all', year = nmy, season = nms, iter = nmi))  
        }   
        else
            if(sl == 'model') slot(a,sl) <- 'PellaTom'
                else
                    if(sl == 'params') { 
                        nmparams <- all.vars(get('PellaTom')()[[2]])
                        nmparams <- nmparams[!(nmparams %in% c('biomass','bprod'))] 
                        slot(a, sl) <- array(dim = c(length(nmparams),ny,ns,ni), dimnames = list(params = nmparams, year = nmy, season = nms, iter = nmi))
                        }
   }

    validObject(a)
    return(a)
}


# BDsim   {{{
# It fills in  the biomass slot for only 1 year and 1 season. 
# Otherwise we would have ny*ns combinations and for simulation we simulate one by one! 
# iter >= 1.

BDsim <- function(object, year = 1, season = 1, iter = 'all')  # year and season either numeric (position) or character (name)
  {

    Dim    <- dim(object@biomass)
    dimnms <- dimnames(object@biomass)
    
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
      
    # Previous season/year.
    if(ss == 1){
        yr0 <- yr-1
        ss0 <- Dim[4]
    }
    else{
        yr0 <- yr
        ss0 <- ss -1
    }    
    
    # model call
    if(length(grep('~', object@model)) == 0)
        model <- eval(call(object@model))[[2]]
    else # character but 'formula' 
        model <- formula(object@model)
    
    model <- as.list(model)[[3]]
     
    # Extract biomass  # numeric[1 OR it]
    datam <- list(biomass = c(object@biomass[,yr0,,ss0,]))

    # Extract catches  # numeric[1 OR it]
    # datam[['catch']] <- c(object@catch[,yr,,ss,])

    # Extract covars.
    for(i in names(object@covar))
        datam[[i]] <-  c(object@covar[[i]][,yr0,, ss0,])
    
    # Extract params
    for(i in dimnames(object@params)[[1]])
        datam[[i]] <-  c(object@params[i,yr0,ss0,])
    
  #  res <- numeric(Dim[6])
    
  #  for(i in 1:Dim[6])
  
  res <- eval(model, datam)
  newB <- object@biomass[,yr0,,ss0,] - object@catch[,yr0,,ss0,] + res*object@uncertainty[,yr0,,ss0,]
    
  object@biomass[,yr,,ss,] <- newB
    
  return(object)

}   # }}}


## Pella-Tomlinson surplus production model
# (from: Jennings et al., 2001. Marine Fisheries Ecology. Blackwell Science, Oxford, UK.)

PellaTom <- function(){

    ## log likelihood, assuming normal log.
    logl <- function( biomass, bprod, r, K, p)
      loglAR1( log(bprod), log(biomass * (r/p) * (1 - (biomass/K)^p)))

    ## model to be fitted
    model <- bprod ~ biomass * (r/p) * (1 - (biomass/K)^p)
    
    # Output
    return( list(logl=logl, model=model) )
}


## iter {{{
setMethod("iter", signature(object="FLBDsim"),
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

