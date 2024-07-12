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

# Weights-at-age based on total biomass and a length frequency distribution

ddwAgeLFD <- function(biol, stknm, year, season, ctrl, covars, ...) {
  
  lfd         <- ctrl[['LFD']]                                                 # Mean length-at-age in a numeric vector.
  a           <- ctrl[['a.lw']]                                                # The "a" parameter of the length-weight relationship.
  dd_mod      <- ctrl[['LFD_model']]                                           # A function that calculates the "b" parameter of the length-weight relationship based on the total biomass.
  excluded.a  <- ctrl[['exc.a']]                                               # By default, the weights of all ages are recalculated, but the weight of specific ages can be fixed by indicating here which one they are.
  
  mx.lfd = matrix(lfd, dim(biol@n)[1], dim(biol@n)[6])
  
  B <- quantSums((biol@wt*biol@n)[,year-1,,season,])[drop=T] 
  
  condb <- predict(dd_mod, data.frame(B = B))                                  # Predict the "b" parameter of the length-weight relationship.
  mx.condb <- matrix(condb, dim(biol@n)[1], dim(biol@n)[6], byrow = TRUE)
  
  wage <- a * mx.lfd^mx.condb/1000 # in tonnes
  
  if(!is.null(excluded.a)){
    excluded.a.pos <- which(biol@range[["min"]]:biol@range[["max"]] %in% excluded.a)  
    wage[excluded.a.pos,] <- biol@wt[excluded.a.pos,year,,season,drop=TRUE]
  }
  
  # wt change
  wt.ref <- biol@wt[,year,,season,drop=TRUE]
  wt.chg <- wage/wt.ref
  
  return(list(wt = wage, wt.chg = wt.chg))
  
}
