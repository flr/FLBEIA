#-------------------------------------------------------------------------------
# FLMetier extension  - Copy FLMetier from FLcore and replace FLMetier by 
# FLMetierExt and FLCatches by FLCatche sExt.
# Dorleta Garcia - 11/08/2010 09:19:28
#-------------------------------------------------------------------------------


## FLMetierExt		{{{
# Check in dimension 4 removed. 
# It produces warnings when quant dimension = 'all'.
validFLMetierExt <- function(object) {
    options(warn = -1)
	# FLQuant slots share dims 1:5 ...
  dnames <- qapply(object, function(x) dimnames(x)[4:5])
	for(i in names(dnames))
		if(!identical(dnames[[i]], dnames[[1]]))
			return('All FLQuant slots must have the same dimensions')

  # ... and are consistent with catches
  catdnames <- lapply(object@catches, function(x)
    qapply(object, function(x) dimnames(x)[4:5]))
  for(i in seq(length=length(catdnames)))
    for(j in names(catdnames[[1]]))
	    if(!identical(catdnames[[i]][[j]], dnames[[1]]))
			  return('All FLQuant slots must have the same dimensions')

  # Year range of FLMetier covers all catches
  catyears <- matrix(unlist(lapply(object@catches, function(x)
    unlist(dims(x)[c('minyear', 'maxyear')]))), byrow=TRUE, ncol=2)
  if(any(dims(object)$minyear < catyears [,1]) |
    any(dims(object)$maxyear > catyears [,2]))
    return('Year range of metier should encompass those of catch(es)')

  # iter is consistent between metier and catches
  if(any(dims(object)$iter != unlist(lapply(object@catches, function(x) dims(x)$iter))))
    return('iter must be 1 or N across all slots and levels')

	return(TRUE)
}

setClass('FLMetierExt',
	representation('FLComp',
		gear='character',
		effshare='FLQuant',
		vcost='FLQuant',
		catches='FLCatchesExt'),
	prototype(name=character(0), desc=character(0),
		range= unlist(list(min=0, max=0, plusgroup=NA, minyear=1, maxyear=1)),
		gear=character(0), catches=new('FLCatchesExt'), effshare=FLQuant(1), vcost=FLQuant(NA)),
	validity=validFLMetierExt)

remove(validFLMetierExt)
# Accesors
# createFLAccesors('FLMetierExt', exclude=c('range', 'catches', 'name', 'desc'))
# }}}

# FLMetierExt
setGeneric('FLMetierExt', function(catches, ...) standardGeneric('FLMetierExt'))

## FLMetier()	{{{
# FLMetier(FLCatch)
setMethod('FLMetierExt', signature(catches='FLCatchExt'),
	function(catches, gear='NA', ...){
	    FLMetierExt(catches=FLCatchesExt(catches), gear=gear, ...)
#	   res <- new('FLMetierExt')
#	   res@catches <- FLCatchesExt(catches)
#	   res@gear    <- gear
	   }
)
# FLMetier(FLCatches)
setMethod('FLMetierExt', signature(catches='FLCatchesExt'),
	function(catches, gear='NA', ...)
    {
		args <- list(...)
    if(length(args) > 0)
    {
      classes <- lapply(args, class)
      # if any in ... is FLQuant
      if(any('FLQuant' %in% classes))
        # take dimnames of first one
        dimn <- dimnames(args[[names(classes[classes %in% 'FLQuant'])[1]]])
    }
    if(!exists('dimn'))
    {
      # generate from FLCatch
      dimn <- dimnames(landings.n(catches[[1]]))
      dimn[[3]] <- 'unique'
      years <- apply(as.data.frame(lapply(catches, function(x) unlist(dims(x)[c(
        'minyear','maxyear')]))), 1, max)
      dimn$year <- as.character(seq(years[1], years[2]))
      dimn[[1]] <- 'all'
    }
    
    # new object
		res <- new('FLMetierExt', catches=catches, gear=gear, effshare=FLQuant(1, dimnames=dimn),
      vcost=FLQuant(NA, dimnames=dimn), range=range(catches))
    # load extra arguments
		if(length(args) > 0)
			for (i in seq(length(args)))
				slot(res, names(args[i])) <- args[[i]]
		return(res)
    }
)
# FLMetier(FLQuant)
setMethod('FLMetierExt', signature(catches='FLQuant'),
	function(catches, gear='NA', ...)
      return(FLMetierExt(FLCatchExt(catches), gear=gear, ...))
)
# FLMetier(missing)
setMethod('FLMetierExt', signature(catches='missing'),
	function(catches, gear='NA', ...)
    FLMetierExt(FLCatchesExt(FLCatchExt(name='NA')), ...)
)	# }}}




# summary	{{{
setMethod('summary', signature(object='FLMetierExt'),
	function(object, ...)
	{
		callNextMethod(object)
		cat("\n")
		cat("Catches: ", "\n")
		for (j in names(object@catches))
			cat("\t", j, ": [", dim(object@catches[[j]]@landings.n),"]\n")
	}
)
# }}}

# trim {{{
setMethod('trim', signature(x='FLMetierExt'),
  function(x, ...)
  {
    x <- callNextMethod()
    x@catches <- lapply(x@catches, trim, ...)
    return(x)
  }
) # }}}

# propagate {{{
setMethod('propagate', signature(object='FLMetierExt'),
  function(object, ...)
  {
    object <- qapply(object, propagate, ...)
    object@catches <- lapply(object@catches, propagate, ...)
    return(object)
  }
) # }}}

## iter {{{
setMethod("iter", signature(obj="FLMetierExt"),
	  function(obj, iter)
	  {
		# FLQuant slots
		names <- names(getSlots(class(obj))[getSlots(class(obj))=="FLQuant"])
		for(s in names) 
		{
			if(dims(slot(obj, s))$iter == 1)
				slot(obj, s) <- iter(slot(obj, s), 1)
			else
				slot(obj, s) <- iter(slot(obj, s), iter)
		}
		# FLCatches
		names <- names(obj@catches)
		for (s in names)
			catches(obj, s) <- iter(catches(obj, s), iter)
		
		return(obj)
	  }
) # }}}

# "[" and "[["             {{{
setMethod("[", signature(x="FLMetierExt", i="ANY", j="missing"),
  function(x, i, drop=FALSE)
  {
	  if (missing(i))
      return(x)
    x@catches <- x@catches[i]
    return(x)
	}
)

setMethod("[[", signature(x="FLMetierExt", i="ANY", j="missing"),
  function(x, i, drop=FALSE)
  {
	  if (missing(i))
      stop("invalid subscript type")
    return(x@catches[[i]])
	}
)  # }}}
