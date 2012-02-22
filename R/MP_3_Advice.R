#-------------------------------------------------------------------------------
#       advice.mp(stocks, fleets.obs, covars, advice,  advice.ctrl, year, season,... )
#       OUTPUT: advice, fleets.obs
#
# Dorleta Garcia
# Created: 20/12/2010 16:21:59
# Changed: 13/01/2011 12:27:18
#-------------------------------------------------------------------------------

advice.mp <- function(stocks, fleets.obs, indices, covars, advice, advice.ctrl, year, season){
   
    stknms <- unique(c(names(stocks), names(indices)))
   
    for(stknm in stknms){
  
  #  if(stknm == 'FAKE') browser()
    
    #cat(stknm,'\n')
    
      #  la siguiente linea era para permitir gestion en diferentes seasons y pluriannual por ahora esta idea estYYY parada.
      #  if(!(year %in% advice.ctrl[[stknm]]$ass.year & season %in% advice.ctrl[[stknm]]$ass.season))  next

     cat('----------------- ', stknm, ' -----------------\n')
     
        advice <- eval(call(advice.ctrl[[stknm]]$HCR, stocks = stocks, covars = covars,  stknm = stknm,
                advice = advice, year = year, indices = indices, advice.ctrl = advice.ctrl))
        
    }
    
    # Apply fleet based advice, this could affect both, the advice itself and the annual quota-share among fleets. 
    # OR overall and fleet especific effort restrictions.
    if(!is.null(advice.ctrl[['fleet']]$HCR))
        advice <- eval(call(advice.ctrl[['fleet']]$HCR, fleet.obs = fleet.obs, 
        stock = stocks, covars = covars, advice = advice, year = year, advice.ctrl = advice.ctrl,...))     
       
    return(advice)

}

fixedAdvice <- function(advice,...){
    return(advice)
}










