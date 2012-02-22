#-------------------------------------------------------------------------------
# FLFleet extension  - Copy FLFleet from FLcore and replace:
#   FLFleet by FLFleetExt, FLMetier(s) by FLMetier(s)Ext and 
#   FLCatch(es) by FLCatch(es) sExt. 
# Dorleta García - 11/08/2010 10:31:07
#-------------------------------------------------------------------------------

## FLFleetExt		{{{
validFLFleetExt <- function(object) {

	# FLQuant slots share dims 3:5 ...
  dnames <- qapply(object, function(x) dimnames(x)[3:5])
	for(i in names(dnames))
		if(!identical(dnames[[i]], dnames[[1]]))
			return('All FLQuant slots must have the same dimensions')

  # ... and are consistent with metiers
  metdnames <- lapply(object@metiers, function(x)
    qapply(object, function(x) dimnames(x)[3:5]))
  for(i in seq(length=length(metdnames)))
    for(j in names(metdnames[[1]]))
	    if(!identical(metdnames[[i]][[j]], dnames[[1]]))
			  return('All FLQuant slots must have the same dimensions')

  # Year range of FLFleetExt covers all metiers
  metyears <- matrix(unlist(lapply(object@metiers, function(x)
    unlist(dims(x)[c('minyear', 'maxyear')]))), byrow=TRUE, ncol=2)

  if(any(dims(object)$minyear < metyears [,1]) |
    any(dims(object)$maxyear > metyears [,2]))
    return('Year range of fleet should encompass those of metier(s)')

  # iter is consistent between fleet and metiers
  if(any(dims(object)$iter != unlist(lapply(object@metiers, function(x) dims(x)$iter))))
    return('iter must be 1 or N across all slots and levels')

  # effshares must add up to one
  #effshs <- lapply(object@metiers, effshare)
  #if(length(effshs) > 1)
  #  for(i in 2:length(effshs))
  #    effshs[[1]] <- effshs[[1]] + effshs[[i]]
  #if(!isTRUE(all.equal(as.vector(effshs[[1]]), rep(1,prod(dim(effshs[[1]]))))))
  #  return('sum of effshare must add up to 1')

	return(TRUE)
}

setClass('FLFleetExt',
	representation('FLComp',
		effort='FLQuant',
		fcost='FLQuant',
		capacity='FLQuant',
		crewshare ="FLQuant",
		metiers='FLMetiersExt'),
	prototype(name=character(0), desc=character(0),
		range= unlist(list(min=0, max=0, plusgroup=NA, minyear=1, maxyear=1)),
		effort=FLQuant(), fcost=FLQuant(), capacity=FLQuant(),
		crewshare=FLQuant(), metiers=FLMetiersExt()),
	validity=validFLFleetExt)
remove(validFLFleetExt)

#invisible(createFLAccesors("FLFleetExt", exclude=c('range', 'effort', 'name', 'desc')))	# }}}

# FLFleetExt
setGeneric('FLFleetExt', function(object, ...) standardGeneric('FLFleetExt'))



# FLFleetExt()		{{{
setMethod('FLFleetExt', signature(object='FLMetiersExt'),
	function(object, ...)
	{
		args <- list(...)
		flqs <- unlist(lapply(args, is, 'FLQuant'))
		if(any(flqs))
			flqs <- FLQuant(NA,
				dimnames=c(dimnames(args[[names(flqs[flqs==TRUE])[1]]])[-6], list(iter=1)))
		else
			flqs <- FLQuant()
		res <- new('FLFleetExt', metiers=object, effort=flqs, fcost=flqs,
			capacity=flqs, crewshare=flqs, range=range(object))

		# extra arguments
		for (i in names(args))
			slot(res, i) <- args[[i]]
		return(res)
	}
)


setMethod('FLFleetExt', signature(object='FLMetierExt'),
	function(object, ...)
	{
		FLFleetExt(FLMetiersExt(met=object), ...)
	}
)
setMethod('FLFleetExt', signature(object='FLCatchesExt'),
	function(object, ...)
	{
		FLFleetExt(FLMetiersExt(FLMetierExt(object)), ...)
	}
)
setMethod('FLFleetExt', signature(object='FLCatchExt'),
	function(object, ...)
	{
		FLFleetExt(FLMetiersExt(FLMetierExt(FLCatchesExt(object))), ...)
	}
)
setMethod('FLFleetExt', signature(object='FLFleetExt'),
	function(object, metier, catch, ...)
	{
    res <- object
    if(!missing(metier))
      res@metiers <- res@metiers[metier]
    if(!missing(catch))

		FLFleetExt(, ...)
	}
)
setMethod('FLFleetExt', signature(object='missing'),
	function(object, ...)
	{
		FLFleetExt(FLMetiersExt(FLMetierExt(FLCatchesExt(FLCatchExt()))), ...)
	}
)	# }}}

# summary	{{{
setMethod('summary', signature(object='FLFleetExt'),
	function(object, ...)
	{
		callNextMethod(object)
		cat("\n")
		cat("Metiers: ", "\n")
		# TODO What happens when object has no metiers/catches? IM 28.08.07
		for (i in names(object@metiers))
		{
			cat("\t", i, ":\n")
			
			for (j in names(object@metiers[[i]]@catches))
				cat("\t\t", j, ": [", dim(object@metiers[[i]]@catches[[j]]@landings.n),"]\n")
		}
	}
)
# }}}

# metier(fl, me)	{{{
setMethod('metier', signature(object='FLFleetExt', metier='ANY'),
	function(object, metier, ...)
		return(object@metiers[[metier]])
)
setReplaceMethod('metier', signature(object='FLFleetExt', metier='ANY', value='FLMetierExt'),
	function(object, metier, ..., value)
	{
		object@metiers[[metier]] <- value
		return(object)
	}
)	# }}}

# FLFleetExt accesors	{{{
createFleetExtAccesors('catch', catch, c(2:5), assigment=FALSE)
createFleetExtAccesors('catch.n', catch.n, c(2:5), assigment=FALSE)
createFleetExtAccesors('catch.wt', catch.wt, c(2:5), assigment=FALSE)
createFleetExtAccesors('catch.sel', catch.sel, c(2:5), assigment=FALSE)
createFleetExtAccesors('catch.q', catch.q)
createFleetExtAccesors('discards', discards)
createFleetExtAccesors('discards.n', discards.n)
createFleetExtAccesors('discards.wt', discards.wt)
createFleetExtAccesors('discards.sel', discards.sel)
createFleetExtAccesors('landings', landings)
createFleetExtAccesors('landings.n', landings.n)
createFleetExtAccesors('landings.wt', landings.wt)
createFleetExtAccesors('landings.sel', landings.sel)
createFleetExtAccesors('price', price)
# }}}

# revenue	{{{
setMethod('revenue', signature('FLCatchExt'),
	function(object)
    if(!all(is.na(landings.n(object))))
      return(quantSums(landings.n(object) * landings.wt(object) * price(object)))
    else
      return(landings(object) * price(object))
)
setMethod('revenue', signature('FLCatchesExt'),
	function(object, catch=unique(names(object)), ...)
		return(lapply(object, revenue))
)
setMethod('revenue', signature('FLMetierExt'),
  function(object, ...)
    return(revenue(object@catches, ...))
)
setMethod('revenue', signature('FLMetiersExt'),
  function(object, metier, catch, ...)
  {
  if(missing(catch) && missing(metier))
    return(TRUE)
  else if(missing(catch))
    revenue(metier(object, metier))
  else if(missing(metier))
    Sums(lapply(object@metiers, revenue))
  else
    return(TRUE)
  }
)
setMethod('revenue', signature('FLFleetExt'),
  function(object, ...)
    return(revenue(object@metiers, ...))
) # }}}

## iter {{{
setMethod("iter", signature(object="FLFleetExt"),
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
		# FLMetiersExt
		names <- names(object@metiers)
		for (s in names)
			metier(object, s) <- iter(metier(object, s), iter)
		
		return(object)
	  }
) # }}}

# catches(fl, me, ca)	{{{
setMethod('catches', signature(object='FLFleetExt'),
	function(object, ...)
		return(catches(object@metiers, ...))
)
setMethod('catches', signature(object='FLMetiersExt'),
	function(object, catch='missing', sum=FALSE, ...)
  {
    # No catch? OK if only one in object
    if(missing(catch))
      if(length(unique(unlist(lapply(object, function(x) names(x@catches))))) == 1)
        catch <- object[[1]]@catches[[1]]@name
      else
        stop('No catch was selected and object holds data for more than one catch')
    
    # identify metiers with this catch.
    idx <- unlist(lapply(object, function(x) any(catchNames(x) == catch)))

    # if index is numeric and only one metier, select from names
    if(length(object) == 1 & is.numeric(catch))
      catch <- catchNames(object)[catch]
    res <- lapply(object[idx], catches, catch=catch)
    
    if(length(res) > 1 && sum==TRUE)
    {
      res[1:2] <- mcf(res[[1]], res[[2]])
      res[[1]] <- addFLCatch(res[[1]], res[[2]])
      if(length(res) > 2)
        for(i in seq(3, length(res)))
        {
          res[[i]] <- mcf(res[[1]], res[[i]])[[2]]
          res[[1]] <- addFLCatch(res[[1]], res[[i]])
        }
      return(FLCatchesExt(res[[1]]))
    }
    return(FLCatchesExt(res))
  }
)
setMethod('catches', signature(object='FLMetierExt'),
	function(object, catch='missing', ...)
  {
		if(missing(catch))
      return(object@catches)
    if (length(catch) == 1)
      return(object@catches[[catch]])
    else
      return(object@catches[catch])
  }
)	# }}}

# catches<-(fl, ca)	{{{
setMethod('catches<-', signature(object='FLMetierExt', value='FLCatchExt'),
	function(object, catch, ..., value)
  {
    object@catches[[catch]] <- value
    return(object)
  }
)
setMethod('catches<-', signature(object='FLMetierExt', value='FLCatchesExt'),
	function(object, catch, ..., value)
  {
    object@catches <- value
    return(object)
  }
) # }}}

# FLMetierExt accesors for FLFleetExt {{{
setMethod('effshare', signature(object='FLMetiersExt'),
  function(object, metier=names(object))
  {
    if(length(metier) == 1)
      return(object[[metier]]@effshare)
    else
      return(FLQuants(lapply(object[metier], effshare)))
  }
)
setMethod('effshare', signature(object='FLFleetExt'),
  function(object, ...)
    return(effshare(object@metiers, ...))
)
setMethod('vcost', signature(object='FLMetiersExt'),
  function(object, metier=names(object))
  {
    if(length(metier) == 1)
      return(object[[metier]]@vcost)
    else
      return(FLQuants(lapply(object[metier], vcost)))
  }
)
setMethod('vcost', signature(object='FLFleetExt'),
  function(object, ...)
    return(vcost(object@metiers, ...))
)
# }}}

## dims {{{
setMethod("dims", signature(obj="FLFleetExt"),
  # Returns a list with different parameters
  function(obj, ...)
	{
		qnames <- names(getSlots(class(obj))[getSlots(class(obj))=="FLQuant"])
		return(list(
      metiers=names(obj@metiers),
      catches=unique(unlist(lapply(obj@metiers, function(x) names(x@catches)))),
      quant = quant(slot(obj, qnames[1])),
      min=min(as.numeric(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) dimnames(x@landings.n)[[1]][1]))))),
      max=max(as.numeric(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) dimnames(x@landings.n)[[1]][dim(x@landings.n)[1]]))))),
      minyear=min(as.numeric(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) dimnames(x@landings.n)[[2]][1]))))),
      maxyear=max(as.numeric(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) dimnames(x@landings.n)[[2]][dim(x@landings.n)[2]]))))),
      unit=unique(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) length(dimnames(x@landings.n)[[3]]))))),
      season=unique(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) length(dimnames(x@landings.n)[[4]]))))),
      area=unique(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) length(dimnames(x@landings.n)[[5]]))))),
      iter=max(unlist(lapply(obj@metiers, function(x) lapply(x@catches,
        function(x) qapply(x, function(x) length(dimnames(x)[[6]]))))))
    ))
    }
)    # }}}

## window    {{{
setMethod("window", signature(x="FLFleetExt"),
	  function(x, start=dims(x)$minyear, end=dims(x)$maxyear, extend=TRUE, frequency=1) {

    # window fleet
    x <- qapply(x, window, start, end, extend, frequency)

    # window metiers
    x@metiers <- lapply(x@metiers, window, start, end, extend, frequency)

    # window catches
    for(i in seq(length(x@metiers)))
      x@metiers[[i]]@catches <- lapply(x@metiers[[i]]@catches, window, start, end, extend, frequency)

		x@range["minyear"] <- start
		x@range["maxyear"] <- end

		return(x)
	}
)	# }}}

## effort		{{{
setMethod("effort", signature(object="FLFleetExt", metier="missing"),
	function(object)
    return(slot(object, "effort")))

setMethod("effort", signature(object="FLFleetExt", metier="character"),
	function(object, metier)
    return(slot(object, "effort") * slot(metier(object, metier), "effshare")))

setReplaceMethod("effort", signature(object="FLFleetExt", value="FLQuant"),
	function(object, value)
  {
		slot(object, "effort") <- value
    return(object)
  })
# }}}

# catchNames {{{
setMethod('catchNames', signature(object='FLCatchesExt'),
  function(object)
  {
    return(unname(unlist(lapply(object, catchNames))))
  }
)
setMethod('catchNames', signature(object='FLMetierExt'),
  function(object)
  {
    return(catchNames(object@catches))
  }
)
setMethod('catchNames', signature(object='FLMetiersExt'),
  function(object)
  {
    return(unique(unlist(lapply(object, catchNames))))
  }
)
setMethod('catchNames', signature(object='FLFleetExt'),
  function(object)
  {
    return(catchNames(object@metiers))
  }
) 
setMethod('catchNames', signature(object='FLFleetsExt'),
  function(object)
  {
    return(unique(unlist(lapply(object, catchNames))))
  }
) # }}}

# trim {{{
setMethod('trim', signature(x='FLFleetExt'),
  function(x, ...)
  {
    x <- callNextMethod()
    x@metiers <- lapply(x@metiers, trim, ...)
    return(x)
  }
) # }}}

# propagate {{{
setMethod('propagate', signature(object='FLFleetExt'),
  function(object, ...)
  {
    object <- qapply(object, propagate, ...)
    object@metiers <- lapply(object@metiers, propagate, ...)
    return(object)
  }
) # }}}

# computeCatch  {{{
setMethod('computeCatch', signature(object='FLCatchExt'),
  function(object)
    return(quantSums(catch.n(object) * catch.wt(object)))
)
setMethod('computeDiscards', signature(object='FLCatchExt'),
  function(object)
    return(quantSums(discards.n(object) * discards.wt(object)))
)
setMethod('computeLandings', signature(object='FLCatchExt'),
  function(object)
    return(quantSums(landings.n(object) * landings.wt(object)))
)

setMethod('computeCatch', signature(object='FLMetierExt'),
  function(object, catch=names(object@catches))
  Sums(lapply(object@catches[catch], computeCatch))
)
setMethod('computeDiscards', signature(object='FLMetierExt'),
  function(object, catch=names(object@catches))
  lapply(object@catches[catch], computeDiscards)
)
setMethod('computeLandings', signature(object='FLMetierExt'),
  function(object, catch=names(object@catches))
  lapply(object@catches[catch], computeLandings)
)

setMethod('computeCatch', signature(object='FLFleetExt'),
  function(object, ...)
  lapply(object@metiers, computeCatch, ...)
)
setMethod('computeDiscards', signature(object='FLFleetExt'),
  function(object, ...)
  lapply(object@metiers, computeDiscards, ...)
)
setMethod('computeLandings', signature(object='FLFleetExt'),
  function(object, ...)
  lapply(object@metiers, computeLandings, ...)
) # }}}

# "[" and "[["             {{{
setMethod("[", signature(x="FLFleetExt", i="ANY", j="missing"),
  function(x, i, drop=FALSE)
  {
	  if (missing(i))
      return(x)
    x@metiers <- x@metiers[i]
    return(x)
	}
)

setMethod("[", signature(x="FLFleetExt", i="ANY", j="ANY"),
  function(x, i, j, drop=FALSE)
  {
    if(!missing(i))
      x <- x[i]
    if(!missing(j))
      x@metiers <- lapply(x@metiers, '[', j)
    return(x)
	}
)

setMethod("[[", signature(x="FLFleetExt", i="ANY", j="missing"),
  function(x, i, drop=FALSE)
  {
	  if (missing(i))
      stop("invalid subscript type")
    return(x@metiers[[i]])
	}
) # }}}

# as.data.frame {{{
setMethod('as.data.frame', signature(x='FLFleetExt', row.names='missing',
  optional='missing'), function(x)
  {
    df <- callNextMethod()
    df <- cbind(df, metier='NA', catch='NA')

    for (i in 1:length(x@metiers))
    {
      df <- rbind(df, cbind(catch='NA', metier=names(x@metiers)[[i]],
        as.data.frame(x@metiers[[i]])))

      for (j in 1:length(x@metiers[[i]]@catches))
      df <- rbind(df, cbind(catch=names(x@metiers[[i]]@catches)[[j]],
        metier=names(x@metiers)[[i]], as.data.frame(x@metiers[[i]]@catches[[j]])))
    }
    return(df)
  }
) # }}}
