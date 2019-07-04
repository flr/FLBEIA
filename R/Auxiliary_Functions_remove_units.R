#-------------------------------------------------------------------------------  
# setUnitsNA function.
# Created: Sonia Sanchez - 2018-11-15
# Changed: 2019-06-25 09:54:55 (ssanchez) - convert function to method
#------------------------------------------------------------------------------- 

# setUnitsNA_method.R - function to remove units from FLQuants of an object
# ./R/setUnitsNA_method.R

# Copyright: AZTI, 2019
# Author: Dorleta Garcia & Sonia Sanchez (AZTI)
#
# Distributed under the terms of the XXX Public Licence (XXX) V.X.X.


# setUnitsNA {{{

#' @title Method setUnitsNA
#' 
#' @description Function to remove the units from the FLQuants of an object
#'
#' @param object An object of class FLBiol, FLBiols, FLFleetExt or FLFleetsExt object.
#' 
#' @return The same object with the units equal to NA in all the FLQuant slots.
#' 
#' @keywords setUnitsNA methods
#' 
#' @examples
#'
#' data(one)
#' 
#' st <- setUnitsNA(oneBio)
#' fl <- setUnitsNA(oneFl)
#' 
#' data(multi)
#' 
#' biols <- setUnitsNA(multiBio)
#' fleets <- setUnitsNA(multiFl)
#' 
#' data(res_flbeia)
#' 
#' res_biols  <- setUnitsNA(multiRes$biols)
#' res_fleets <- setUnitsNA(multiRes$fleets)
#' 
#

#----------------------------------------------------------------------------

#' @rdname setUnitsNA
#' @aliases setUnitsNA,FLBiol-method
# setUnitsNA(object='FLBiol',...) {{{
setMethod('setUnitsNA', signature(object='FLBiol'),
          function(object)
          {
            
            for(sl in c('n', 'wt', 'spwn', 'm')){ 
              units(slot(object, sl))[] <- NA
              
            }
            
            for(sl in c('fec', 'mat')){
              units(slot(object, sl)[[sl]])[] <- NA
            }
            
            return(object)
          }
)
# }}}

#' @rdname setUnitsNA
#' @aliases setUnitsNA,FLBiols-method
# setUnitsNA(object='FLBiols',...) {{{
setMethod('setUnitsNA', signature(object='FLBiols'),
          function(object)
          {
            
            for (st in names(object))
              object[[st]] <- setUnitsNA(object[[st]])
            
            return(object)
            
          }
)
# }}}

#' @rdname setUnitsNA
#' @aliases setUnitsNA,FLFleetExt-method
# setUnitsNA(object='FLFleetExt',...) {{{
setMethod('setUnitsNA', signature(object='FLFleetExt'),
          function(object)
          {
            
            for(sl in c("effort", "fcost", "capacity", "crewshare")){ 
              units(slot(object, sl))[] <- NA
            }
            
            for(mt in names(object@metiers)){
              for(sl in c('effshare', 'vcost')){
                units(slot(object@metiers[[mt]], sl))[] <- NA
              }
              for(ct in names(object@metiers[[mt]]@catches)){
                for(sl in c("landings", "landings.n", "landings.wt", "landings.sel", 
                            "discards", "discards.n", "discards.wt", "discards.sel", 
                            "catch.q",  "price")){
                  units(slot(object@metiers[[mt]]@catches[[ct]], sl))[] <- NA
                }
              }
            }
            
            return(object)
          }
)
# }}}

#' @rdname setUnitsNA
#' @aliases setUnitsNA,FLFleetsExt-method
# setUnitsNA(object='FLFleetsExt',...) {{{
setMethod('setUnitsNA', signature(object='FLFleetsExt'),
          function(object)
          {
            
            for (fl in names(object))
              object[[fl]] <- setUnitsNA(object[[fl]])
            
            return(object)
            
          }
)
# }}}

