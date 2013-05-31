#-------------------------------------------------------------------------------
#       advice.mp(stocks, fleets.obs, covars, advice,  advice.ctrl, year, season,... )
#       OUTPUT: advice, fleets.obs
#
# Dorleta Garcia
# Created: 20/12/2010 16:21:59
# Changed: 13/01/2011 12:27:18
# Changes: 2012-06-15 13:26:13  Sonia Sánchez - for allowing assessment in different seasons and multiannual advice
#-------------------------------------------------------------------------------

advice.mp <- function(stocks, fleets.obs, indices, covars, advice, advice.ctrl, year, season, stknm){
   
   cat('----------------- ', stknm, ' -----------------\n')
   
      advice <- eval(call(advice.ctrl[[stknm]]$HCR.model, stocks = stocks, covars = covars,  stknm = stknm,
              advice = advice, year = year, indices = indices, advice.ctrl = advice.ctrl))
      
    # Apply fleet based advice, this could affect both, the advice itself and the annual quota-share among fleets. 
    # OR overall and fleet especific effort restrictions.
    if(!is.null(advice.ctrl[['fleet']]$HCR.model))
        advice <- eval(call(advice.ctrl[['fleet']]$HCR.model, fleets.obs = fleets.obs, 
        stock = stocks, covars = covars, advice = advice, year = year, advice.ctrl = advice.ctrl,...))     
       
    return(advice)

}

fixedAdvice <- function(advice,...){
    return(advice)
}










