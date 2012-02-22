## createFLeetAccesors  {{{
createFleetExtAccesors <- function(slot, fun, level=c(1:5), assigment=TRUE, class='FLQuant')
  {
	# replacement function
  if(assigment == TRUE)
  {
	# FLCatch
  if(1 %in% level)
	eval(substitute(setReplaceMethod(SLOT, signature(object='FLCatchExt', value=class),
		function(object, value) {
			slot(object, SLOT) <- value
			return(object)}), list(SLOT=slot)))
	# FLCatches
  if(2 %in% level)
	eval(substitute(setReplaceMethod(SLOT, signature(object='FLCatchesExt', value=class),
		function(object, catch, value) {
			slot(object[[catch]], SLOT) <- value
			return(object)
		}),list(SLOT=slot)))
	# FLMetier
  if(3 %in% level)
	eval(substitute(setReplaceMethod(SLOT, signature(object='FLMetierExt', value=class),
		function(object, catch, value) {
			slot(object@catches[[catch]], SLOT) <- value
			return(object)
		}),list(SLOT=slot)))
	# FLMetiers
  if(4 %in% level)
	eval(substitute(setReplaceMethod(SLOT, signature(object='FLMetiersExt', value=class),
		function(object, metier, catch, value) {
			slot(object[[metier]]@catches[[catch]], SLOT) <- value
			return(object)
		}), list(SLOT=slot)))
	# FLFleet
  if(5 %in% level)
	eval(substitute(setReplaceMethod(SLOT, signature(object='FLFleetExt', value=class),
		function(object, metier, catch, value) {
			slot(object@metiers[[metier]]@catches[[catch]], SLOT) <- value
			return(object)
		}), list(SLOT=slot)))
  }

	# accesor functions
	# FLCatch
  if(1 %in% level)
	eval(substitute(setMethod(SLOT, signature(object='FLCatchExt'),
		function(object)
			return(slot(object, SLOT))), list(SLOT=slot)))
	# FLCatches
  if(2 %in% level)
	eval(substitute(setMethod(SLOT, signature(object='FLCatchesExt'),
		function(object, catch='missing') {
			if(missing(catch))
				return(lapply(object, SLOT))
			else
				return(FUN(object[[catch]]))}),list(SLOT=slot, FUN=fun)))
	# FLMetier
  if(3 %in% level)
	eval(substitute(setMethod(SLOT, signature(object='FLMetierExt'),
		function(object, ...)
				return(FUN(object@catches, ...))), list(SLOT=slot, FUN=fun)))
	# FLMetiers
  if(4 %in% level)
	eval(substitute(setMethod(SLOT, signature(object='FLMetiersExt'),
		function(object, metier='missing', catch='missing', ...) {
      # nothing
			if (missing(metier) && missing(catch))
				stop('Either metier or catch must be specified')
      # metier
			else if(!missing(metier) && missing(catch))
				return(FUN(object[[metier]], ...))
      # catch
			else if(missing(metier) && !missing(catch))
      {
				res <- FLQuants()
				for(i in names(object))
        {
          if (catch %in% names(object[[i]]@catches))
  					res[[i]] <- FUN(object[[i]], catch=catch, ...)
        }
				return(res)
      # both
			} else
				return(FUN(object[[metier]], catch=catch, ...))}), list(SLOT=slot, FUN=fun)))
	# FLFleet
  if(5 %in% level)
	eval(substitute(setMethod(SLOT, signature(object='FLFleetExt'),
		function(object, ...)
				return(FUN(object@metiers, ...))), list(SLOT=slot, FUN=fun)))
}   # }}}
