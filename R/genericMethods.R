#-------------------------------------------------------------------------------  
# genericMethods.
# Changed: 
#------------------------------------------------------------------------------- 

# genericMethods - S4 generics
# FLCore/R/genericMethods.R

# Copyright: AZTI, 2019
# Author: Sonia Sanchez (AZTI) <ssanchez@azti.es>)
#
# Distributed under the terms of the XXX Licence V.x.x.
# 
# Maintainer: Sonia Sanchez, AZTI



# # globalVariables(c("qname"))
# # 
# # -- OVERLOADED methods/functions {{{
# 
# setGeneric("AIC", useAsDefault = stats::AIC)
# 
# # }}}



# -- CONSTRUCTORS, documented with each class {{{


# # XXX
# #' @rdname XXX
# #' @aliases XXX XXX-methods
# setGeneric("XXX", function(object, ...) standardGeneric("XXX"))
# 
# # }}}


# # -- ACCESSORS {{{
# 
# #' @rdname XXX-class
# #' @aliases XX XX-methods
# setGeneric("XX",function(object, ...) standardGeneric("XX"))
# 
# # }}}



# -- METHODS


# setUnitsNA {{{
#' @rdname setUnitsNA
#' @aliases setUnitsNA setUnitsNA-methods
setGeneric("setUnitsNA", function(object)
  standardGeneric("setUnitsNA")) 

# }}}



