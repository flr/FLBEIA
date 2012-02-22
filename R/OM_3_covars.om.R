#-------------------------------------------------------------------------------
#       covars.om(covars, biols, fleets, covars.ctrl, year, season,... )
# 
#   Output: Updated FLBiols 
# 
#
# Dorleta Garcia
# Created: 17/01/2011 15:55:21
# Changed: 17/01/2011 15:55:26
#-------------------------------------------------------------------------------

covars.om <- function(biols, fleets, covars, advice, covars.ctrl, year, season){


    cvnms <- names(covars)
    
    for(cv in cvnms){
        # population dynamic model
        dyn.model <- covars.ctrl[[cv]]$dyn.model   
        
        res <- eval(call(dyn.model, covars = covars, biols = biols, fleets = fleets, advice = advice, year = year, season = season, ctrl = covars.ctrl, cvnm = cv))
        
        covars       <- res[['covars']]
        fleets       <- res[['fleets']]
        biols        <- res[['biols']]
    }
    
    return(list(covars = covars, fleets = fleets, biols = biols))

}


#-------------------------------------------------------------------------------
# NoDynCovar - No dynamics, return the object without modifying.
#-------------------------------------------------------------------------------

fixedCovar <-  function(covars, biols, fleets, cvnm,...){

    covars[[cvnm]] <- covars[[cvnm]]
    
    return(list(covars = covars, fleets = fleets, biols = biols))
}









