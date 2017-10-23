#-------------------------------------------------------------------------------
#                   STOCK RECRUITMENT FUNCTIONS
#
# 20/09/2011 12:08:50
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMARK: '...' in the arguments of the functions are neccesary in order to be
#   generalistic inside 'biols.om' function ('eval(call(...)').
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


hockstick <- function () {
  
    logl <- function(a, b, rec, ssb) {
        loglAR1(log(rec), FLQuant( log(ifelse(c(ssb) <= b, a * c(ssb), a * b)), 
                                   dimnames = dimnames(ssb)))
    }
  
    model <- rec ~ ifelse(c(ssb) <= b, a * c(ssb), a * b)
  
    initial <- structure(function(rec, ssb) {
        return(FLPar(a = median(c(rec/ssb), na.rm = TRUE), b = median(c(ssb), na.rm = TRUE)))
    }, lower = rep(0, 0), upper = rep(Inf, 2))
  
  return(list(logl = logl, model = model, initial = initial))
}



# Redfish Recruitment model, idea by Benjamin Planque coding by Dorleta Garcia
redfishRec <- function(recPrevY,sigma, minrec, maxrec){
   
    
    it <- length(sigma)
    rec <- numeric(it)
    
    for(i in 1:it){
         k <- 0
        while(k<1){ 
            rec[i] <- recPrevY[i] + rnorm(n=1,mean=0,sd=sigma[i])
            if((rec[i] >= minrec[i]) & (rec[i] <= maxrec[i])) k <- k+1
        }
    }
    return(rec)
}

redfishRecModel <- function(){
    logl <- NA
    model <- rec ~ ifelse(ssb < alpha, ssb/alpha, 1)*redfishRec(rec.prevY,sigma, minrec, maxrec)
    initial <- NA
    return(list(logl = logl, model = model, initial = initial))
}

# Ricker models, whith a covariate in the exponent:
# - ANCHOVY with predation from SARDINE
aneRec_pil <- function () 
{
  
  logl <- function(a, b, c, rec, ssb, ssb.pil) 
    loglAR1(log(rec), log(a * ssb * exp(-b * ssb + c * ssb.pil)))
  
  model <- rec ~ a * ssb * exp(-b * ssb + c * ssb.pil)
  
  return(list(logl = logl, model = model))
  
}
# - SARDINE with predation from ANCHOVY
pilRec_ane <- function () 
{
  
  logl <- function(a, b, c, rec, ssb, ssb.ane) 
    loglAR1(log(rec), log(a * ssb * exp(-b * ssb + c * ssb.ane)))
  
  model <- rec ~ a * ssb * exp(-b * ssb + c * ssb.ane)
  
  return(list(logl = logl, model = model))
  
}

# - constant recruitment
constRec <- function (a, rec, ssb) 
{
  logl <- NA
  model <- rec ~ a
  return(list(logl = logl, model = model))  
}


# Ricker model, whith a cyclic term:
ricker.cyclic <- function () 
{
  
  logl <- function(a, b, c, d, p, rec, ssb, year) 
    loglAR1(log(rec), log(a * ssb * exp(-b * ssb + c * sin(2 * pi * year/p) + d * cos(2 * pi * year/p))))
  
  model <- rec ~ a * ssb * exp(-b * ssb +  c * sin(2 * pi * year/p) + d * cos(2 * pi * year/p))
  
  return(list(logl = logl, model = model))
  
}

# Ricker with years (yr.low=T) with low recruitment (rec.low)
ricker.low <- function(){
  logl <- NA
  model <- rec ~ yr.low * rec.low + (1-yr.low) * (a * ssb * exp(-b * ssb)) # i.e. Ricker 
  initial <- NA
  return(list(logl = logl, model = model, initial = initial))
}

