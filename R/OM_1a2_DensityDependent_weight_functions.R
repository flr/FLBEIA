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
# (specific for Bay of Biscay anchovy with 3 age classes - 0 to 3+)

ddwAgeCa <- function(biol, stknm, year, season, ctrl, covars, ...) {
  
  require(nleqslv)
  
  pars  <- ctrl$params[,,year,season,] # array[npar,nage,nyr,ns,nit]
  covnm <- ctrl$covnm                  # name of the covariate
  
  wt.ref <- biol@wt[,year,,season,drop=TRUE]
  ssna   <- (n(biol) * mat(biol) * exp(-spwn(biol) * m(biol)))[,year,,season,drop=TRUE] # exclude age0 (inmature)
  sst    <- covars[[covnm]][stknm,year,,season,drop=TRUE]
  
  # ssb.ref <- quantSums(ssna * wt.ref[-1])
  
  wfun <- function(x) {
    
    # x[i] is log(wage_i)
    # x[i] = a_i + b_i * log(ssb) + c_i * sst
    #      = a_i + b_i * log(sum_j[ssna_j * exp(x[j])]) + c_i * sst 
    #         with ssna_j = @n_j * @mat_j * exp(-@m_j * @spwn_j)
    
    y <- numeric(3)
    
    logssb <- log(ssna[1] * exp(x[1]) + ssna[2] * exp(x[2]) +
                    ssna[3] * exp(x[3]) + ssna[4] * exp(x[4]))
    
    y[1] <- x[1] - pars["b",1] * logssb - pars["a",1] - pars["c",1] * sst
    y[2] <- x[2] - pars["b",2] * logssb - pars["a",2] - pars["c",2] * sst
    y[3] <- x[3] - pars["b",3] * logssb - pars["a",3] - pars["c",3] * sst
    y[4] <- x[4] - pars["b",4] * logssb - pars["a",4] - pars["c",4] * sst
    
    return(y)
  }
  
  lwest <- nleqslv(log(wt.ref), fn = wfun, control = list(btol=1e-08, delta="newton"))
  
  if (lwest$termcd != 1)
    print(paste0(stknm," weights calculation message: ",lwest3$message))
  
  if (any(wfun(lwest$x)>1e-08))
    stop(paste("Issues with mean weight calculation for", stknm))
  
  wage <- exp(lwest$x)
  
  # check
  logssb <- sum(ssna * wage)
  di <- numeric(4)
  for (i in 1:4) {
    di[i] <- wage[i] - exp(pars["a",i] + pars["b",i] * logssb + pars["c",i] * sst)
  }
  if (any(round(di-wage,8) != 0))
    stop(paste("Issues with density-dependent weights-at-age for",stknm))
  
  # wt change
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

