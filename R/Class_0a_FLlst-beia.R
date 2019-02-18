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

#' Class FLCatchesExt
#' 
#' A list of \code{FLCatchExt} objects.
#'
#' @name FLCatchesExt
#' @aliases FLCatchesExt-class 
#' @docType class
#' @section Slots: \describe{
#'   \item{.Data}{Internal S4 data representation, of class \code{list}.}
#'   \item{desc}{As textual description of the object contents}
#'   \item{lock}{Can the object be extended/trimmed? \code{TRUE} or \code{FALSE}.}
#'   \item{names}{A character vector for the element names} }
# @template FLlst-constructors
#' @param object An object of class FLCatchExt, list or missing.
#' @param ... Other objects to be assigned by name to the class slots.
#' 
#' @author The FLBEIA Team
#' @seealso \code{\link{FLlst}}, \code{\link[base]{list}},
#'   \code{\link[base]{vector}}, \code{\link{FLCatchExt}}
#' @keywords classes
#'


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
#' @rdname FLCatchesExt
setClass("FLCatchesExt", contains="FLlst", 
	validity=vFLCs
)

# constructor
setGeneric("FLCatchesExt", function(object, ...){
	standardGeneric("FLCatchesExt")
	}
)

#' @aliases FLCatchesExt,ANY
#' @rdname FLCatchesExt
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

#' @aliases FLCatchesExt,missing
#' @rdname FLCatchesExt
setMethod("FLCatchesExt", "missing", function(...){
	if(missing(...)){
		new("FLCatchesExt")
	} else { 
		lst <- list(...)
		new("FLCatchesExt", lst)
	}
})

#' @aliases FLCatchesExt,list
#' @rdname FLCatchesExt
setMethod("FLCatchesExt", "list", function(object){
	new("FLCatchesExt", object)
})


#' @aliases is.FLCatchesExt
#' @rdname FLCatchesExt
setGeneric("is.FLCatchesExt", function(object, ...){
	standardGeneric("is.FLCatchesExt")
	}
)

#' @aliases is.FLCatchesExt,ANY
#' @rdname FLCatchesExt
setMethod("is.FLCatchesExt", "ANY", function(object, ...){
	identical(is(object)[1],"FLCatchesExt")
})

#! FLMetiersExt

#' Class FLMetiersExt
#' 
#' A list of \code{FLMetierExt} objects.
#'
#' @name FLMetiersExt
#' @aliases FLMetiersExt-class 
#' @docType class
#' @section Slots: \describe{
#'   \item{.Data}{Internal S4 data representation, of class \code{list}.}
#'   \item{desc}{As textual description of the object contents}
#'   \item{lock}{Can the object be extended/trimmed? \code{TRUE} or \code{FALSE}.}
#'   \item{names}{A character vector for the element names} }
# @template FLlst-constructors
#' @param object An object of class FLMetierExt, list or missing.
#' @param ... Other objects to be assigned by name to the class slots.
#' 
#' @author The FLBEIA Team
#' @seealso \code{\link{FLlst}}, \code{\link[base]{list}},
#'   \code{\link[base]{vector}}, \code{\link{FLMetierExt}}
#' @keywords classes
#'


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
#' @rdname FLMetiersExt
setClass("FLMetiersExt", contains="FLlst", 
	validity=vFLMs
)

# constructor
setGeneric("FLMetiersExt", function(object, ...){
	standardGeneric("FLMetiersExt")
	}
)

#' @aliases FLMetiersExt,ANY
#' @rdname FLMetiersExt
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
#' @rdname FLMetiersExt
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
#' @rdname FLMetiersExt
setMethod("FLMetiersExt", "list", function(object){
	new("FLMetiersExt", object)
})

# is
#' @aliases FLMetiersExt,is.FLMetiersExt
#' @rdname FLMetiersExt
setGeneric("is.FLMetiersExt", function(object, ...){
	standardGeneric("is.FLMetiersExt")
	}
)

#' @aliases FLMetiersExt,is.FLMetiersExt,ANY
#' @rdname FLMetiersExt
setMethod("is.FLMetiersExt", "ANY", function(object, ...){
	identical(is(object)[1],"FLMetiersExt")
})

#! FLFleetsExt

#' Class FLFleetsExt
#' 
#' A list of \code{FLFleetExt} objects.
#'
#' @name FLFleetsExt
#' @aliases FLFleetsExt-class 
#' @docType class
#' @section Slots: \describe{
#'   \item{.Data}{Internal S4 data representation, of class \code{list}.}
#'   \item{desc}{As textual description of the object contents}
#'   \item{lock}{Can the object be extended/trimmed? \code{TRUE} or \code{FALSE}.}
#'   \item{names}{A character vector for the element names} }
# @template FLlst-constructors
#' @param object An object of class FLFleetExt, list or missing.
#' @param ... Other objects to be assigned by name to the class slots.
#' 
#' @author The FLBEIA Team
#' @seealso \code{\link{FLlst}}, \code{\link[base]{list}},
#'   \code{\link[base]{vector}}, \code{\link{FLFleetExt}}
#' @keywords classes
#'


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
#' @rdname FLFleetsExt
setClass("FLFleetsExt", contains="FLlst",
	validity=vFLFs
)

# constructor
setGeneric("FLFleetsExt", function(object, ...){
	standardGeneric("FLFleetsExt")
	}
)

#' @aliases FLFleetsExt,ANY
#' @rdname FLFleetsExt
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
#' @rdname FLFleetsExt
setMethod("FLFleetsExt", "missing", function(...){
	if(missing(...)){
		new("FLFleetsExt")
	} else { 
		lst <- list(...)
		new("FLFleetsExt", lst)
	}
})

#' @aliases FLFleetsExt,list
#' @rdname FLFleetsExt
setMethod("FLFleetsExt", "list", function(object){
	new("FLFleetsExt", object)
})

# is
#' @aliases is.FLFleetsExt
#' @rdname FLFleetsExt
setGeneric("is.FLFleetsExt", function(object, ...){
	standardGeneric("is.FLFleetsExt")
	}
)

#' @aliases is.FLFleetsExt,ANY
#' @rdname FLFleetsExt
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
