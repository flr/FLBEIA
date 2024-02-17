#-------------------------------------------------------------------------------
#                   DENSITY-DEPENDENT WEIGHT FUNCTIONS
#
# 14/02/2024
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMARK: '...' in the arguments of the functions are necessary in order to be
#   generalist inside 'ASPG_DDW' function ('eval(call(...)').
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Input variables from ASPG_DWW call:
# biol = biol, stknm = stknm, year = year, season = season, ctrl = ddw.ctrl, covars = covars


#-------------------------------------------------------------------------------
# ddwAgeCa(biol, SR, fleets, biol.control)
# - OUTPUT: list(wt = wage, wt.chg = wt.chg) - vector with estimated weight at age values and relative change
#-------------------------------------------------------------------------------

# Weights-at-age based on SSB and a covariate

ddwAgeCa <- function(biol, stknm, year, season, ctrl, covars, ...) {
  
  pars  <- ctrl$params[,,year,season,] # array[npar,nage,nyr,ns,nit]
  covnm <- ctrl$covnm                  # name of the covariate
  
  logssb  <- log(ssb(biol)[,year,,season,drop=TRUE])
  wage    <- exp(pars["a",] + pars["b",] * logssb + pars["c",] * covars[[covnm]][stknm,year,,season,drop=TRUE])
  
  # wt change
  wt.ref <- biol@wt[,year,,season,]
  wt.chg <- wage/wt.ref
  
  return(list(wt = wage, wt.chg = wt.chg))
  
}


#-------------------------------------------------------------------------------
# ddwAgeLFD(biol, SR, fleets, biol.control)
# - OUTPUT: list(wt = wage, wt.chg = wt.chg) - vector with estimated weight at age values and relative change
#-------------------------------------------------------------------------------

# Weights-at-age based on SSB (linear model for estimating LW b parameter) and A LFD

ddwAgeLFD <- function(biol, stknm, year, season, ctrl, covars, ...) {
  
  lfd       <- ctrl[['LFD']]
  a         <- ctrl[['a.lw']]
  lbins     <- as.numeric(colNames(lfd))
  
  B <- quantSums((biol@wt*biol@n)[,year-1])[drop=T] #! DG needs to consider season dimension
  
  condF <- predict(LW_lm, data.frame(biomass = B))  #! DG requires: biols.ctrl[[stknm]][['ddw.ctrl']][['LW_lm']]
  
  wy <- a*(lbins)^condF
  
  wt. <- rowSums(sweep(lfd, 2, wy, "*"))  #! DG needs to consider also season dimension
  wt  <- biol@wt[,year,,season,]
  
  return(list(wt = wt., wt.chg = wt./wt))
  
}

