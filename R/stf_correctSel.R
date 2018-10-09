# stf - prepare an object for a Short Term Forecast
# FLAssess/R/stf.R

# Copyright 2003-2008 FLR Team. Distributed under the GPL 2 or later

## stf(FLStock) {{{

stf_correctSel<-  function(object, nyears=3, wts.nyears=3, fbar.nyears=wts.nyears, f.rescale=FALSE,
    arith.mean=TRUE, na.rm=TRUE, end=dims(object)$maxyear + nyears)
  {
    dims <- dims(object)

    # check nyears and end match
    if(missing(nyears))
      nyears <- as.numeric(end) - dims$maxyear
    else if(dims$maxyear + nyears != end)
      stop("'nyears' and 'end' do not match: ", dims$maxyear + nyears, " vs. ", end)
    years      <- ac((dims$maxyear+1):end)
    wts.years  <- ac(seq(dims$maxyear-wts.nyears+1, dims$maxyear))
    fbar.ages  <- ac(range(object, 'minfbar'):range(object, 'maxfbar'))
    fbar.years <- ac(seq(dims$maxyear-fbar.nyears+1, dims$maxyear))
    years.catch_G0<- which(object@catch> 1e-2)
    fbar.years <- tail(years.catch_G0,3)
    nit <- dim(object@catch)[6]
    fbar.years.iters <- matrix(0,nrow=3,ncol=nit)
    
    for(ii in 1:nit){
      years.catch_G0<- which(iter(object@catch,ii)> 1e-5)
      fbar.years.iters[,ii] <- tail(years.catch_G0,3)
    }
    

    # arith or geometric
    if(arith.mean)
      fmean <- mean
    else  
      fmean <- function(x) exp(mean(log(x)))

    # window object
    res <- window(object, end=end)

    # average slots
    # *.wt, mat, m and *.spwn as average over wts.years
    for (i in c('catch.wt', 'landings.wt', 'discards.wt', 'stock.wt', 'mat', 'm', 'harvest.spwn', 'm.spwn')){
      flq<- apply(slot(res, i)[,wts.years], c(1,3:6),fmean, na.rm=na.rm)
      for (j in years)
         slot(res, i)[,j] <-flq
      }

    # landings.n and discards.n as proportions of wts.years
    for (i in years)
       slot(res, 'discards.n')[,i] <- apply(slot(res, 'discards.n')[, wts.years]/slot(res, 'catch.n')[, wts.years], c(1,3:6), mean)
    slot(res, 'landings.n')[,years] <- 1 - slot(res, 'discards.n')[,years]

    # first f get the dimensions and after harvest as mean over fbar.years.iter
    f <-apply(slot(res, 'harvest')[, wts.years], c(1,3:6), fmean, na.rm=na.rm)
    for(ii in 1:nit){
      iter(f,ii)  <-apply(iter(slot(res, 'harvest'),ii)[,fbar.years.iters[,ii]], c(1,3:6), fmean, na.rm=na.rm)}
      
    for (i in years)
       slot(res, 'harvest')[,i] <- f

    # f.rescale
    if(f.rescale == TRUE)
    {
      # mean f over fbar ages and years
      fbar <- yearMeans(apply(slot(res, 'harvest')[fbar.ages, wts.years], c(2:6), mean,
        na.rm=na.rm))
      lastfbar <- apply(slot(res, 'harvest')[fbar.ages, tail(wts.years,1)], 3:6, mean,
                        na.rm=na.rm)
      
      for(ii in 1:nit){
        iter(fbar,ii)[]  <- mean(apply(iter(slot(res, 'harvest'),ii)[fbar.ages, fbar.years.iters[,ii]], c(2:6), mean,
                                     na.rm=na.rm))
        iter(lastfbar,ii)[] <- apply(iter(slot(res, 'harvest'),ii)[fbar.ages, tail(fbar.years.iters[,ii],1)], 3:6, mean,
                          na.rm=na.rm)
        iter(slot(res, 'harvest'),ii)[, years] <- sweep(iter(slot(res, 'harvest'),ii)[, years], 3:6,iter(fbar,ii), '/')
        iter(slot(res, 'harvest'),ii)[, years] <- sweep(iter(slot(res, 'harvest'),ii)[, years], 3:6, iter(lastfbar,ii), '*')
        }
      

      # # fbar for last REAL year
      # lastfbar <- apply(slot(res, 'harvest')[fbar.ages, tail(fbar.years,1)], 3:6, mean,
      #   na.rm=na.rm)

      # divide by fbar and multiply by lastfbar
      # slot(res, 'harvest')[, years] <- sweep(slot(res, 'harvest')[, years], 3:6, fbar, '/')
      # slot(res, 'harvest')[, years] <- sweep(slot(res, 'harvest')[, years], 3:6, lastfbar, '*')
    }
    return(res)
  }
 # }}}

