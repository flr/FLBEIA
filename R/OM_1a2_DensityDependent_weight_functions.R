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
  
  na <- dim(biol@wt)[1]
  it <- dim(biol@wt)[6]
  
  pars  <- ctrl$params[,,year,season,] # array[npar,nage,1,1,nit]
  covnm <- ctrl$covnm                  # name of the covariate
  
  wt.ref <- biol@wt[,year,,season,drop=TRUE] # [na,nit] or [na]
  ssna   <- (n(biol) * mat(biol) * exp(-spwn(biol) * m(biol)))[,year,,season,drop=TRUE] # [na,nit] or [na] # to exclude age0 (immature)
  sst    <- covars[[covnm]][stknm,year,,season,drop=TRUE]  # [nit] 
  
  wage <- wt.ref * NA
  
  # ssb.ref <- quantSums(ssna * wt.ref[-1])
  
  wfun <- function(x, param, ssna) {
    
    # x[i] is log(wage_i)
    # x[i] = a_i + b_i * log(ssb) + c_i * sst
    #      = a_i + b_i * log(sum_j[ssna_j * exp(x[j])]) + c_i * sst 
    #         with ssna_j = @n_j * @mat_j * exp(-@m_j * @spwn_j)
    
    y <- numeric(3)
    
    logssb <- log(ssna[1] * exp(x[1]) + ssna[2] * exp(x[2]) +
                    ssna[3] * exp(x[3]) + ssna[4] * exp(x[4]))
    
    y[1] <- x[1] - param["b",1] * logssb - param["a",1] - param["c",1] * sst
    y[2] <- x[2] - param["b",2] * logssb - param["a",2] - param["c",2] * sst
    y[3] <- x[3] - param["b",3] * logssb - param["a",3] - param["c",3] * sst
    y[4] <- x[4] - param["b",4] * logssb - param["a",4] - param["c",4] * sst
    
    return(y)
  }
  
  if (it == 1) {
    
    lwest <- nleqslv(log(wt.ref), fn = wfun, param=pars, ssna=ssna, control = list(btol=1e-08, delta="newton"))
    
    if (lwest$termcd != 1)
      print(paste0(stknm," weights calculation message: ",lwest3$message))
    
    if (any(wfun(lwest$x, param=pars, ssna=ssna)>1e-08))
      stop(paste("Issues with mean weight calculation for", stknm))
    
    wage <- exp(lwest$x)
    
    # # check
    # logssb <- sum(ssna * wage)
    # di <- numeric(na)
    # for (a in 1:na)
    #   di[a] <- wage[a] - exp(pars["a",a] + pars["b",a] * logssb + pars["c",a] * sst)
    # if (any(round(di-wage,4) != 0))
    #   stop(paste("Issues with density-dependent weights-at-age for",stknm))
    
  } else {
    
    for (i in 1:it) {
      
      lwest <- nleqslv(log(wt.ref[,i]), fn = wfun, param=pars[,,i], ssna=ssna[,i], control = list(btol=1e-08, delta="newton"))
      
      if (lwest$termcd != 1)
        print(paste0(stknm," weights calculation message: ",lwest3$message))
      
      if (any(wfun(lwest$x, param=pars[,,i], ssna=ssna[,i])>1e-08))
        stop(paste("Issues with mean weight calculation for", stknm))
      
      wage[,i] <- exp(lwest$x)
    }
    
    # # check
    # logssb <- apply(ssna * wage, 2, sum)
    # di <- array(NA, dim=c(na,it))
    # for (a in 1:na)
    #   di[a,] <- wage[a,] - exp(pars["a",a,] + pars["b",a,] * logssb + pars["c",a,] * sst)
    # if (any(round(di-wage,4) != 0))
    #   stop(paste("Issues with density-dependent weights-at-age for",stknm))
    
  }
    
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
