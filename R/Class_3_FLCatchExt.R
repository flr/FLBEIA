#-------------------------------------------------------------------------------
# FLCatchExt - Extension of original FLCatch object, it includes 2 new slots:
# alpha and beta: Coob-Douglas function parameters.
#
# Dorleta Garcia - Azti Tecnalia - 21/07/2010
#   Change: 06/06/2011 (dga) extend FLCatch to make flexible quant dimension in alpha, beta and catch.q.
#-------------------------------------------------------------------------------

## FLCatchExt               {{{
validFLCatchExt <- function(object)
{
    names <- names(getSlots('FLCatchExt')[getSlots('FLCatchExt')=="FLQuant"])
    nits  <- sort(unique(unlist(qapply(object, function(x) dims(x)$iter))))

  if (length(nits)>2) return(paste("All FLQuant must either have same number of iters or '1 & n'"))

	for(i in names)
	{
		# all dimnames but iter, age and unit are the same
		if(!identical(unlist(dimnames(object@landings.n)[c(2,4:5)]),
			unlist(dimnames(slot(object, i))[c(2,4:5)])))
			return(paste('All elements must share dimensions 2,4 and 5: Error in FLCatchExt', i))
	}
	for (i in names[!names%in%c('landings', 'discards', 'catch.q', 'alpha', 'beta')])
	{
		# quant is n
		if(!identical(unlist(dimnames(object@landings.n)[1]),
			unlist(dimnames(slot(object, i))[1])))
			return(paste('All elements must share quant names: Error in FLCatchExt', i))
	}

	for (i in c('alpha', 'beta', 'catch.q'))
	{
		# quant is n
		if(!identical(unlist(dimnames(object@alpha)), unlist(dimnames(slot(object, i)))))
			return('"alpha", "beta" and "catch.q" must share dimensions: Error in FLCatchExt')
	}
	# alpha, beta and catch.q have the same dimensions, check that alpha has correct dimension in quant.
#	if(dim(object@alpha)[3] > 1 & dim(object@alpha)[3] != dim(object@alpha)[4]) return('Wrong unit dimension for "alpha", "beta" and "catch.q"')
	if(dim(object@alpha)[1] > 1 & dim(object@alpha)[1] != dim(object@landings.n)[1]) return('Wrong quant dimension for "alpha", "beta" and "catch.q"')

	for (i in c('landings', 'discards'))
	{
		# quant is 1
		if(any(dim(slot(object, i))[c(1,3,5)] != 1))
			return(paste('Wrong dimensions for slot ', i, 'in FLCatchExt'))
	}
	return(TRUE)
}

#' 
#' @name FLCatchExt
#' @aliases FLCatchExt-class FLCatchExt FLCatchExt-methods
#' FLQuant,FLCatchExt-method
#' 
#' @title  FLCatchExt class and the methods to construct it.
#'
#' @description It extends the FLCatch class defined in FLFleet package. 
#'    The FLCatch class includes two extra slots alpha and beta used in the Cobb-Douglas production functions.
#' 
#' @details The FLCatchExt object contains a representation of the catch of a fish stock as constructed for the purposes of fleet dynamic modelling. 
#'    This includes information on removals (i.e. landings and discards), selectivity, weight, price and catch production parameters (catchability and elasticities).
#'    
#' @param object,x An object of class FLQuant, missing or FLCatchExt.
#' @param range Numerical vector with min, max, plusgroup, minyear and maxyear elements as in FLStock object.
#' @param name The name of the stock.
#' @param desc The description of the object.
#' @param i,j,k,l,m,n subindices
#' @param drop logical. Should the dimesions be dropped?
#' @param ... Other objects to be assigned by name to the class slots.
#' @param e1,e2 FLCatchExt objects, where e2 is incorporated into e1 (see addFLCatch).
#' @param value Value or values to be assigned to the particular FLQuant or FLCatchExt slot.
#' 
#' @slot landings An FLQuant with the total landings in weight of the stock.
#' @slot landings.n An FLQuant with the landings in numbers at age of the stock.
#' @slot landings.wt An FLQuant with the weight at age of the landings.
# @slot landings.sel An FLQuant with the retention ogive of the metier for this stock. 
#                    Elements must be between 0 and 1. landings.sel = 1-discards.sel. 
#' @slot discards An FLQuant with the total discards in weight of the stock.
#' @slot discards.n An FLQuant with the discards in numbers at age of the stock.
#' @slot discards.wt An FLQuant with the weight at age of the discards.
#' @slot landings.sel,discards.sel An FLQuant with the landing/discard ogive of the metier for this stock 
#'                                 (i.e. landings.sel corresponds to the proportion of catches landed).
#'                                 Elements must be between 0 and 1, with discards.sel = 1-landings.sel. 
#' @slot catch.q An FLQuant with the catchability at age of the stock for the corresponding metier. 
#'               This is the catchability used in the catch production model.
#' @slot price An FLQuant with the price at age of the stock.
#' @slot alpha An FLQuant with the elasticity parameter at age of the stock for the corresponding metier. 
#'             This is one of the parameters used in the catch production model.
#' @slot beta An FLQuant with the elasticity parameter at age of the stock for the corresponding metier. 
#'            This is one of the parameters used in the catch production model.
#' @slot name The name of the stock.
#' @slot desc A description of the object.
#' @slot range The range as in other FLR objects: c("min","max","plusgroup","minyear","maxyear").
#'
#' 
#' @return The constructors return an object of class FLCatchExt.
#'  
#' 
#' 
#' 
setClass("FLCatchExt",
    representation(
		'FLComp',
            landings    = "FLQuant", landings.n   = "FLQuant",
            landings.wt = "FLQuant", landings.sel = "FLQuant",
            discards    = "FLQuant", discards.n   = "FLQuant",
            discards.wt = "FLQuant", discards.sel = "FLQuant",
            catch.q     = "FLQuant", price        = "FLQuant",
            alpha       = "FLQuant", beta         = "FLQuant"),
    prototype=prototype(
		name		= character(0),
		desc		= character(0),
	  range       = as.numeric(c(min=NA, max=NA, plusgroup=NA,
			minyear=NA, maxyear=NA)),
    landings = new("FLQuant"), landings.n = new("FLQuant"),
    landings.wt = new("FLQuant"), landings.sel = new("FLQuant"),
    discards = new("FLQuant"), discards.n  = new("FLQuant"),
    discards.wt = new("FLQuant"), discards.sel= new("FLQuant"),
    catch.q     = new("FLQuant"), price = new("FLQuant"),
    alpha = new("FLQuant"), beta = new("FLQuant")),
	validity=validFLCatchExt
)
remove(validFLCatchExt) # }}}


# FLCatchExt()                {{{
setGeneric('FLCatchExt', function(object, ...)
		standardGeneric('FLCatchExt'))

# TODO Fix size of input objects and validity
#' @aliases FLCatchExt,FLQuant-method
#' @rdname FLCatchExt
setMethod('FLCatchExt', signature(object='FLQuant'),
	function(object, range='missing', name='NA', desc=character(0), ...) {
		# initial objects
		flq <- FLQuant(NA, dimnames=dimnames(object))
		flqa <- quantSums(flq)
		flqau <- unitSums(quantSums(flq))
		dims <- dims(flq)
		args <- list(...)

    # construct range
		if(missing(range))
			range <- c(min=dims$min, max=dims$max, plusgroup=NA,
				minyear=dims$minyear, maxyear=dims$maxyear)

		# output object
		res <- new('FLCatchExt', range=range, name=name, desc=desc,
			landings.n=flq, landings.wt=flq, landings.sel=flq, landings=flqau,
			discards.n=flq, discards.wt=flq, discards.sel=flq, discards=flqau,
			catch.q=flq, price=flq, alpha = flq, beta = flq)
		# Load given slots
		for(i in names(args))
			slot(res, i) <- args[[i]]
		return(res)
	}
)

#' @aliases FLCatchExt-missing
#' @rdname FLCatchExt
setMethod('FLCatchExt', signature(object='missing'),
	function(...)
  {
		# get arguments & select first full FLQuant
		args <- list(...)
		args <- args[lapply(args, class) == 'FLQuant']

		flqs <- args[!(names(args) %in% c('landings', 'discards'))]

		# select full flquant, or small flquant, or create dimnames
		if(length(flqs) > 0)
			dimnames <- dimnames(flqs[[1]])
		else if(length(args) > 0)
			dimnames <- dimnames(args[[1]])
		else
			dimnames <- dimnames(FLQuant())
		return(FLCatchExt(FLQuant(dimnames=dimnames), ...))
	}
)	# }}}


## computeLandings	{{{
# @rdname FLCatchExt
setMethod("computeLandings", signature(object="FLCatchExt"),
	function(object, na.rm=TRUE) {
        res <- quantSums(landings.n(object)*landings.wt(object), na.rm=na.rm)
        units(res) <- paste(units(landings.n(object)), units(landings.wt(object)))
        return(res)
 	}
)	# }}}

## computeDiscards	{{{
setMethod("computeDiscards", signature(object="FLCatchExt"),
	function(object, na.rm=TRUE) {
        res <- quantSums(discards.n(object)*discards.wt(object), na.rm=na.rm)
        units(res) <- paste(units(discards.n(object)), units(discards.wt(object)))
        return(res)
 	}
)	# }}}

# '['       {{{
#' @rdname FLCatchExt
#' @aliases [,FLCatchExt,ANY,ANY,ANY-method
setMethod('[', signature(x='FLCatchExt'),
	function(x, i, j, k, l, m, n, ..., drop=FALSE)
  {
    dn <- dimnames(landings.n(x))

		if (missing(i))
			i  <-  seq(1, length(dn[1][[1]]))
		if (missing(j))
			j  <-  dn[2][[1]]
   	if (missing(k))
   		k  <-  dn[3][[1]]
		if (missing(l))
			l  <-  dn[4][[1]]
		if (missing(m))
			m  <-  dn[5][[1]]
		if (missing(n))
			n  <-  dn[6][[1]]

    # catch.q  dim[1] = 1 => dim[4] = 1
    if(dim(slot(x, 'catch.q'))[1] == 1)
      slot(x, 'catch.q') <- slot(x, 'catch.q')[1,j,k,1,m,n, drop=FALSE]
    else
      slot(x, 'catch.q') <- slot(x, 'catch.q')[i,j,k,1,m,n, drop=FALSE]

    # full quants
	  quants <- list("landings.n", "landings.wt", "landings.sel",
      "discards.n","discards.wt","discards.sel","price")
    for(q in quants)
      slot(x, q) <- slot(x, q)[i,j,k,l,m,n, drop=FALSE]

    # no-quant quants
	  quants <- list("landings", "discards", "alpha", "beta")
    for(q in quants)
      slot(x, q) <- slot(x, q)[1,j,k,1,m,n, drop=FALSE]

    # range
    x@range['min'] <- as.numeric(dn[[1]][1])
    x@range['max'] <- as.numeric(dn[[1]][length(dn[[1]])])
    x@range['minyear'] <- as.numeric(dn[[2]][1])
    x@range['maxyear'] <- as.numeric(dn[[2]][length(dn[[2]])])

    return(x)
  }
)   # }}}

## "[<-"            {{{
#' @rdname FLCatchExt
#' @aliases [<-,FLCatchExt,ANY,ANY,ANY-method
setMethod("[<-", signature(x="FLCatchExt", value="FLCatchExt"),
	function(x, i, j, k, l, m, n, ..., value)
	{
    dn <- dimnames(landings.n(x))

		if (missing(i))
			i  <-  dn[1][[1]]
		if (missing(j))
			j  <-  dn[2][[1]]
   	if (missing(k))
   			k  <-  dn[3][[1]]
		if (missing(l))
			l  <-  dn[4][[1]]
		if (missing(m))
			m  <-  dn[5][[1]]
		if (missing(n))
			n  <-  dn[6][[1]]

	  quants <- list("catch.q", "landings.n", "landings.wt", "landings.sel",
      "discards.n","discards.wt","discards.sel","price")
    for(q in quants)
      slot(x, q)[i,j,k,l,m,n] <- slot(value, q)

    quants <- list("landings", "discards", "alpha", "beta")
    for(q in quants)
      slot(x, q)[1,j,k,l,m,n] <- slot(value,q)

 		return(x)
	}
)   # }}}


# addFLCatch
setGeneric('addFLCatch', function(e1, e2, ...)
		standardGeneric('addFLCatch'))

# catchNames
setGeneric('catchNames', function(object, ...)
		standardGeneric('catchNames'))
#' @aliases catchNames<-
#' @rdname FLCatchExt 
setGeneric('catchNames<-', function(object, ..., value)
		standardGeneric('catchNames<-'))
		
# addFLCatch for FLCatch {{{
#' @aliases addFLCatch,FLCatchExt,FLCatchExt FLCatchExt-methods
#' @rdname FLCatchExt 
 setMethod('addFLCatch', signature(e1='FLCatchExt', e2='FLCatchExt'),
  function(e1, e2)
  {
    # add
    qnames <- c('landings', 'landings.n', 'discards', 'discards.n')
    for(i in qnames)
      slot(e1, i) <- slot(e1, i) + slot(e2, i)
    # mean weighted by catch
    qnames <- c('landings.sel', 'discards.sel')
    for(i in qnames)
      slot(e1, i) <- slot(e1, i) + slot(e2, i)
    # mean weighted by effshare
    return(e1)
  }
) # }}}

# setPlusGroup  {{{
setMethod('setPlusGroup', signature(x='FLCatchExt', plusgroup='numeric'),
	function(x, plusgroup, na.rm=FALSE)
	{
	#check plusgroup valid
	if (!missing(plusgroup))
     x@range["plusgroup"]<-plusgroup
  if(x@range["plusgroup"] > x@range["max"])
		 return("Error : plus group greater than oldest age")

  pg.range <- as.character(x@range["max"]:x@range["plusgroup"])

	for (i in c("landings.wt", "landings.sel", "price"))
     slot(x,i)[as.character(x@range["plusgroup"])]<-quantSums(slot(x,i)[pg.range]*x@landings.n[pg.range])/quantSums(x@landings.n[pg.range])

	for (i in c("discards.wt", "discards.sel"))
     slot(x,i)[as.character(x@range["plusgroup"])]<-quantSums(slot(x,i)[pg.range]*x@discards.n[pg.range])/quantSums(x@discards.n[pg.range])

  x@landings.n[as.character(x@range["plusgroup"])]<-quantSums(x@landings.n[pg.range])
  x@discards.n[as.character(x@range["plusgroup"])]<-quantSums(x@discards.n[pg.range])

  x<-x[as.character(x@range["min"]:x@range["plusgroup"])]

  x@range["max"]<-x@range["plusgroup"]

	return(x)
	}
)# }}}

# catchNames {{{
#' @aliases catchNames,FLCatchExt
#' @rdname FLCatchExt 
setMethod('catchNames', signature(object='FLCatchExt'),
  function(object)
    return(object@name))

#' @aliases catchNames<-,FLCatchExt,character
#' @rdname FLCatchExt 
setReplaceMethod('catchNames', signature(object='FLCatchExt', value='character'),
  function(object, value)
  {
    object@name <- value
    return(object)
  }
) # }}}

## trim     {{{
setMethod("trim", signature("FLCatchExt"), function(x, ...){

	args <- list(...)

  rng<-range(x)


  c1 <- args[[quant(x@landings.n)]]
	c2 <- args[["year"]]
	c3 <- args[["unit"]]
	c4 <- args[["season"]]
	c5 <- args[["area"]]
	c6 <- args[["iter"]]

  # FLQuants with quant
  names <- getSlotNamesClass(x, 'FLQuant')
  for (name in names)
  {
    if(name %in% c('landings', 'discards'))
      slot(x,name) <- trim(slot(x,name), year=c2, unit=c3, season=c4,
        area=c5, iter=c6)
        
    if(name %in% c('alpha', 'beta'))
      slot(x,name) <- trim(slot(x,name), year=c2, season=c4,
        area=c5, iter=c6)
    # catch.q
    else if(name == 'catch.q')
    {
      if(dim(slot(x, 'catch.q'))[1] == 1)
        slot(x, 'catch.q') <- trim(slot(x, 'catch.q'), year=c2,
          season=c4, area=c5, iter=c6)
      else
        slot(x, 'catch.q') <- trim(slot(x, 'catch.q'), ...)
    }
    else
    {
      slot(x,name) <- trim(slot(x,name), ...)
    }
  }

  # range
  if (length(c1) > 0) {
    x@range["min"] <- c1[1]
    x@range["max"] <- c1[length(c1)]
    if (rng["max"] != x@range["max"])
        x@range["plusgroup"] <- NA
  }
  if (length(c2)>0 ) {
    x@range["minyear"] <- c2[1]
    x@range["maxyear"] <- c2[length(c2)]
  }
	return(x)

}) # }}}

# catch et al {{{
# catch
setMethod('catch', signature(object='FLCatchExt'),
  function(object)
  {
    res <- landings(object) + discards(object)
    if (units(discards(object)) == units(landings(object)))
		  units(res) <- units(discards(object))
    else
      warning("units of discards and landings do not match")
    return(res)
  }
)

# catch.n
setMethod('catch.n', signature(object='FLCatchExt'),
  function(object)
  {
    res <- landings.n(object) + discards.n(object)
    if (units(discards.n(object)) == units(landings.n(object)))
		  units(res) <- units(discards.n(object))
    else
      warning("units of discards.n and landings.n do not match")
    return(res)
  }
)

# catch.wt
setMethod('catch.wt', signature(object='FLCatchExt'),
  function(object, method='n')
  {
    if(method == 'n')
    {
      idx <- landings.n(object) + discards.n(object)
      idx[idx == 0]  <- 1

      res <- (landings.wt(object) * landings.n(object) +
          discards.wt(object) * discards.n(object)) / idx
    } else if (method == 'sel')
    {
      idx <- landings.sel(object) + discards.sel(object)
      idx[idx == 0]  <- 1
      res <- (landings.wt(object) * landings.sel(object) +
          discards.wt(object) * discards.sel(object)) / idx
    }
    if (units(discards.wt(object)) == units(landings.wt(object)))
				units(res) <- units(discards.wt(object))
    return(res)
  }
)

# catch.sel
setMethod('catch.sel', signature(object='FLCatchExt'),
  function(object)
  {
    return(landings.sel(object) + discards.sel(object))
  }
) # }}}
