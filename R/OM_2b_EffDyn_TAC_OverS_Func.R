#-------------------------------------------------------------------------------
#                   TAC OVERSHOOT FUNCTIONS
#
#   This functions are intended to be used within effort dynamic functions.
#
# INPUT: fleets, TAC (numeric, the TAC for the season and fleet, i.e, the quota share),  fleets.ctrl
# OUTPUT: numeric. (The percent of overshoot).
# 
#
# Dorleta Garcia
# Created: 16/03/2012 13:42:54
#-------------------------------------------------------------------------------

TAC.OS.triangCond <- function(fleets, TAC, fleets.ctrl, flnm, stknm, year, season){
    
    Cmax  <- c(apply(apply(catchWStock.f(fleets[[flnm]],stknm)[,1:(year-1),,season], c(2,6), sum), 6, max, na.rm=T)) # [it]
    TACst <- TAC[stknm,]
    it <- length(TACst)
    triang.pars <-  fleets.ctrl[[flnm]][[stknm]][['TAC.OS.triangCond.params']] 
    alpha <- rtriangle(it,triang.pars[['min']], triang.pars[['max']], triang.pars[['mode']])
    alpha <- ifelse(TACst*alpha > Cmax, Cmax/TACst, alpha)

    return(alpha)
}



