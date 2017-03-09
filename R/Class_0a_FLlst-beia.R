#=====================================================================
#
#! Class: FLlst & FL*s
#! Date: 13/04/2007
#! Version: 0.1-0
#
# \name{ myclass-class }
# \docType{class}
# \alias{ myclass-class }
# \title{  - class}
# \description{
#	Represents this using that. 
# }
# \section{Slots}{
#	\describe{
#		\item{\code{slot1}:}{Object of class \code{"???"} ...}
#		\item{\code{slot2}:}{Object of class \code{"???"} ...}
# 	}
# }
# \section{Methods}{
# Type \code{showMethods(classes="myclass")} at the R prompt for a complete list of methods which are available for this class.
# Useful methods include:
#	\describe{
#		\item{\code{mymeth1}:}{ blahblahblah }
#		\item{\code{mymeth2}:}{ blahblahblah }
#	}
# }
# \seealso{
#	\code{\link{ myclass2-class}}
# }
# \section{Creating Objects from the Class}{
#	new("myclass")
# }
# \note{
#	If relevant ...	
# }
# \author{
#	me 
#	\email{me@myself.org}
#	\url{www.myself.org}
# }
# \keyword{classes}
# \keyword{classes}
#
#! ToDo:
#
#! References (bibtex):
#
#! Notes:
#
#=====================================================================

# getPlural               
getPlural <- function(object)
{
  switch(class(object),
    'FLCatchExt'='FLCatchesExt',
    'FLMetierExt'='FLMetiersExt',
    'FLFleetExt'='FLFleetsExt',
    'list'
    )
}


#! FLCatchesExt

# validity
vFLCs <- function(object){
	# Make sure the list contains all items of the same class
	for(i in 1:length(object)){
		if(!is(object[[i]], "FLCatchExt")) stop("Components must be FLCatchExt")	
	}
	# Everything is fine
	return(TRUE)
}

# class
#' @aliases FLCatchesExt
#' @rdname FLCatchExt
setClass("FLCatchesExt", contains="FLlst", 
	validity=vFLCs
)

# constructor
setGeneric("FLCatchesExt", function(object, ...){
	standardGeneric("FLCatchesExt")
	}
)

#' @aliases FLCatchesExt, ANY
#' @rdname FLCatchExt
setMethod("FLCatchesExt", signature(object="ANY"), function(object, ...){
	lst1 <- list(...)

	nlst <- length(lst1)
	lst <- list()
	length(lst) <- nlst + 1
	lst[[1]] <- object
	lst[-1] <- lst1

	# IM 20.08.07 get names
	names <- c(object@name, unlist(lapply(lst1, function(x) x@name)))
	attr(lst, "names") <- names
	attr(lst, "lock") <- TRUE
	new("FLCatchesExt", lst)
})

#' @aliases FLCatchesExt, missing
#' @rdname FLCatchExt
setMethod("FLCatchesExt", "missing", function(...){
	if(missing(...)){
		new("FLCatchesExt")
	} else { 
		lst <- list(...)
		new("FLCatchesExt", lst)
	}
})

#' @aliases FLCatchesExt, list
#' @rdname FLCatchExt
setMethod("FLCatchesExt", "list", function(object){
	new("FLCatchesExt", object)
})


#' @aliases FLCatchesExt, is.FLCatchesExt
#' @rdname FLCatchExt
setGeneric("is.FLCatchesExt", function(object, ...){
	standardGeneric("is.FLCatchesExt")
	}
)

#' @aliases FLCatchesExt, is.FLCatchesExt,ANY
#' @rdname FLCatchExt
setMethod("is.FLCatchesExt", "ANY", function(object, ...){
	identical(is(object)[1],"FLCatchesExt")
})

#! FLMetiersExt

# validity
vFLMs <- function(object){
	# Make sure the list contains all items of the same class
	for(i in 1:length(object)){
		if(!is(object[[i]], "FLMetierExt")) stop("Components must be FLMetierExt")	
	}
	# Everything is fine
	return(TRUE)
}

# class
#' @aliases FLMetiersExt
#' @rdname FLMetierExt
setClass("FLMetiersExt", contains="FLlst", 
	validity=vFLMs
)

# constructor
setGeneric("FLMetiersExt", function(object, ...){
	standardGeneric("FLMetiersExt")
	}
)

#' @aliases FLMetiersExt,ANY
#' @rdname FLMetierExt
setMethod("FLMetiersExt", signature(object="ANY"), function(object, ...){
	lst1 <- list(...)
	nlst <- length(lst1)
	lst <- list()
	length(lst) <- nlst + 1
	lst[[1]] <- object
	lst[-1] <- lst1

  if(is.null(names(lst1)))
    names(lst1) <- rep("", nlst)
  names(lst) <- c(object@gear, names(lst1))
  if(any(names(lst) == ""))
  names(lst)[names(lst)==""] <- unlist(lapply(lst[names(lst)==""], function(x) x@gear))

	attr(lst, "lock") <- TRUE
	new("FLMetiersExt", lst)
})

#' @aliases FLMetiersExt,missing
#' @rdname FLMetierExt
setMethod("FLMetiersExt", "missing", function(...){
	if(missing(...)){
		new("FLMetiersExt")
	} else { 
		lst <- list(...)
    if(any(names(lst) == ""))
      names(lst)[names(lst)==""] <- unlist(lapply(lst[names(lst)==""], function(x) x@gear))
		new("FLMetiersExt", lst)
	}
})

#' @aliases FLMetiersExt,list
#' @rdname FLMetierExt
setMethod("FLMetiersExt", "list", function(object){
	new("FLMetiersExt", object)
})

# is
#' @aliases FLMetiersExt,is.FLMetiersExt
#' @rdname FLMetierExt
setGeneric("is.FLMetiersExt", function(object, ...){
	standardGeneric("is.FLMetiersExt")
	}
)

#' @aliases FLMetiersExt,is.FLMetiersExt,ANY
#' @rdname FLMetierExt
setMethod("is.FLMetiersExt", "ANY", function(object, ...){
	identical(is(object)[1],"FLMetiersExt")
})

#! FLFleetsExt

# validity
vFLFs <- function(object){
	# Make sure the list contains all items of the same class
	for(i in 1:length(object)){
		if(!is(object[[i]], "FLFleetExt")) stop("Components must be FLFleetExt")	
	}
	# Everything is fine
	return(TRUE)
}

# class
#' @aliases FLFleetsExt
#' @rdname FLFleetExt
setClass("FLFleetsExt", contains="FLlst",
	validity=vFLFs
)

# constructor
setGeneric("FLFleetsExt", function(object, ...){
	standardGeneric("FLFleetsExt")
	}
)

#' @aliases FLFleetsExt,ANY
#' @rdname FLFleetExt
setMethod("FLFleetsExt", signature(object="ANY"), function(object, ...){
	lst1 <- list(...)
	nlst <- length(lst1)
	lst <- list()
	length(lst) <- nlst + 1
	lst[[1]] <- object
	lst[-1] <- lst1
	new("FLFleetsExt", lst)
})

#' @aliases FLFleetsExt,FLFleetExt-missing
#' @rdname FLFleetExt
setMethod("FLFleetsExt", "missing", function(...){
	if(missing(...)){
		new("FLFleetsExt")
	} else { 
		lst <- list(...)
		new("FLFleetsExt", lst)
	}
})

#' @aliases FLFleetsExt,list
#' @rdname FLFleetExt
setMethod("FLFleetsExt", "list", function(object){
	new("FLFleetsExt", object)
})

# is
#' @aliases is.FLFleetsExt
#' @rdname FLFleetExt
setGeneric("is.FLFleetsExt", function(object, ...){
	standardGeneric("is.FLFleetsExt")
	}
)

#' @aliases is.FLFleetsExt,ANY
#' @rdname FLFleetExt
setMethod("is.FLFleetsExt", "ANY", function(object, ...){
	identical(is(object)[1],"FLFleetsExt")
})


#! FLBDSRsims

# validity

# class
# @rdname FLBDSRsim
# setClass("FLBDSRsims", contains="FLlst", 
# 	validity=vFLBDSRsims
# )
# 
# # constructor
# setGeneric("FLBDSRsims", function(object, ...){
# 	standardGeneric("FLBDSRsims")
# 	}
# )

# # @rdname FLBDSRsim
# setMethod("FLBDSRsims", signature(object="ANY"), function(object, ...){
# 	lst1 <- list(...)
# 	nlst <- length(lst1)
# 	lst <- list()
# 	length(lst) <- nlst + 1
# 	lst[[1]] <- object
# 	lst[-1] <- lst1
# 
#   if(is.null(names(lst1)))
#     names(lst1) <- rep("", nlst)
#   names(lst) <- c(object@name, names(lst1))
#   if(any(names(lst) == ""))
#   names(lst)[names(lst)==""] <- unlist(lapply(lst[names(lst)==""], function(x) x@name))
# 
# 	attr(lst, "lock") <- TRUE
# 	new("FLBDSRsims", lst)
# })
# 
# # @rdname FLBDSRsim
# setMethod("FLBDSRsims", "missing", function(...){
# 	if(missing(...)){
# 		new("FLBDSRsims")
# 	} else { 
# 		lst <- list(...)
#     if(any(names(lst) == ""))
#       names(lst)[names(lst)==""] <- unlist(lapply(lst[names(lst)==""], function(x) x@name))
# 		new("FLBDSRsims", lst)
# 	}
# })
# 
# # @rdname FLBDSRsim
# setMethod("FLBDSRsims", "list", function(object){
# 	new("FLBDSRsims", object)
# })
# 
# # is
# # @rdname FLBDSRsim
# setGeneric("is.FLBDSRsims", function(object, ...){
# 	standardGeneric("is.FLBDSRsims")
# 	}
# )
# 
# # @rdname FLBDSRsim
# setMethod("is.FLBDSRsims", "ANY", function(object, ...){
# 	identical(is(object)[1],"FLBDSRsims")
# })
