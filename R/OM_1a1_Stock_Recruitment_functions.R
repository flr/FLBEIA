#-------------------------------------------------------------------------------
#                   STOCK RECRUITMENT FUNCTIONS
#
# 20/09/2011 12:08:50
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMARK: '...' in the arguments of the functions are neccesary in order to be
#   generalistic inside 'biols.om' function ('eval(call(...)').
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Segreg <- function ()
{
    logl <- function(a, b, rec, ssb) {
        loglAR1(log(rec), FLQuant(log(ifelse(c(ssb) <= b, a *
            c(ssb), a * b)), dimnames = dimnames(ssb)))
    }
    model <- rec ~ ifelse(c(ssb) <= b, a * c(ssb), a *
        b)
    initial <- structure(function(rec, ssb) {
        return(FLPar(a = median(c(rec/ssb)), b = median(c(ssb))))
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