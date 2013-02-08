#-------------------------------------------------------------------------------
#       covars.om(covars, biols, fleets, covars.ctrl, year, season,... )
# 
#   Output: Updated FLBiols 
# 
#
# Dorleta Garcia
# Created: 17/01/2011 15:55:21
# Changed: 17/01/2011 15:55:26
# Changed: 2012-07-24 17:39:04 - Sonia Sánchez (AZTI)
#               Added: - variable SRs in covars.om input and output
#                      - ssb.get function
#-------------------------------------------------------------------------------

covars.om <- function(biols, SRs, fleets, covars, advice, covars.ctrl, year, season){


    cvnms <- names(covars)
    
    for(cv in cvnms){
        # population dynamic model
        dyn.model <- covars.ctrl[[cv]]$dyn.model   
        
        res <- eval(call(dyn.model, covars = covars, biols = biols, SRs = SRs, fleets = fleets, advice = advice, year = year, season = season, ctrl = covars.ctrl, cvnm = cv))
        
        covars       <- res[['covars']]
        fleets       <- res[['fleets']]
        biols        <- res[['biols']]
        SRs          <- res[['SRs']]
    }
    
    return(list(covars = covars, fleets = fleets, biols = biols, SRs = SRs))

}


#-------------------------------------------------------------------------------
# NoDynCovar - No dynamics, return the object without modifying.
#-------------------------------------------------------------------------------

fixedCovar <-  function(covars, biols, SRs, fleets, cvnm,...){

    covars[[cvnm]] <- covars[[cvnm]]
    
    return(list(covars = covars, fleets = fleets, biols = biols, SRs = SRs))
}


#-------------------------------------------------------------------------------
# ssb.get - gets a stock ssb in the specified season
#-------------------------------------------------------------------------------

ssb.get <-  function( covars, biols, SRs, fleets, year, season, ctrl, cvnm,...){
  
  stknm   <- ctrl[[cvnm]]$ssb.stock
  ss      <- ctrl[[cvnm]]$spwn.sson
  sr.st   <- ctrl[[cvnm]]$sr.covar
  sr.spwn <- which(SRs[[sr.st]]@proportion[,year,,,,1]==1)
  
  if ( !stknm %in% names(biols))
    stop( "Covariable '", cvnm, "' not available for '", stknm, "' stock")
    
  if ( !cvnm %in% names(SRs[[sr.st]]@covar))
    stop( "Covariable '", cvnm, "' not defined in the '", sr.st, "' stock-recrutiment relationship")
    
  if (ss == season) {
    covars[[cvnm]][,year,,season,] <- ssb(biols[[stknm]])[,year,,season,]
    SRs[[sr.st]]@covar[[cvnm]][,year,,sr.spwn,] <- covars[[cvnm]][,year,,season,]
  }
  
  return(list(covars = covars, fleets = fleets, biols = biols, SRs = SRs))
}

