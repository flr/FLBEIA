#----------------------------------------------------------------------------
#'  Function to remove the units from the FLQuants of FLBiol or FLFleetExt objects 
#'
#' @param obj FLBiol or FLFleetExt object
#' 
#' @return The same object with the units equal to NA in all the FLQuant slots.
#' 
#' 
#' @examples
#'\dontrun{
#'
#' data(one)
#' 
#' fl <- setUnitsNA(oneBio)
#' fl <- setUnitsNA(oneFl)
#'}
# Dorleta Garcia
# 2018/11/15
#----------------------------------------------------------------------------


setUnitsNA <- function(obj){
  if(class(obj) == 'FLBiol'){
    for(sl in c('n', 'wt', 'spwn', 'm')){ 
      units(slot(obj, sl))[] <- NA
      
    }
    for(sl in c('fec', 'mat')){
      units(slot(obj, sl))[[sl]] <- NA
    }
  }
  
  if(class(obj) == 'FLFleetExt'){
    for(sl in c("effort", "fcost", "capacity", "crewshare")){ 
      units(slot(obj, sl))[] <- NA
    }
    for(mt in names(obj@metiers)){
      for(sl in c('effshare', 'vcost')){
        units(slot(obj@metiers[[mt]], sl))[] <- NA
      }
      for(ct in names(obj@metiers[[mt]]@catches)){
        for(sl in c("landings", "landings.n", "landings.wt", "landings.sel", 
                    "discards", "discards.n", "discards.wt", "discards.sel", 
                    "catch.q",  "price")){
          units(slot(obj@metiers[[mt]]@catches[[ct]], sl))[] <- NA
      }
    }
  }}
  
  return(obj)
}