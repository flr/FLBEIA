###################################################33
#
#               SCA 
#
# Developed in JRC meeting 1/02/2017
#..................................................


sca2flbeia <- function(stock, indices, control=control,covars=covars){
  
  # Remove the variability, otherwise sca takes de 'uncertainty' as a weight.
  indices <- FLIndices(lapply(indices, function(x) {
                                    x@index.var[] <- NA
                                    x <- trim(x, age = 0:8)
                                    return(x)}))
  
  stock@landings.n[stock@landings.n == 0] <- 0.1
  stock@catch.n <- stock@landings.n + stock@discards.n
    
  stock <- setPlusGroup(stock, 9)
  
  if(control$test==TRUE){
    fit0 <- sca(stock, indices)
    }else{
    if(!"fmod" %in% names(control)){
    print("Using default fmod sca settings")
    fmod <- getFmod(stock, dfm=c(2/3, 2/3))
    }else{
    fmod <- control$fmod
    }
    if(!"qmod" %in% names(control)){
      print("Using default qmod sca settings")
      qmod <- getQmod(indices)
    }else{
      qmod <- control$qmod
    }
    if(!"srmod" %in% names(control)){
      print("Using default srmod sca settings")
      fit0 <- sca(stock, indices, fmodel=fmod, qmodel=qmod)
    }else{
      srmod <- control$srmod
      fit0 <- sca(stock, indices, fmodel=fmod, qmodel=qmod,srmodel=srmod)
    }
  }
  
  ## convergence diagnostic (quick and dirty)
  # maxgrad <- fit0@fitSumm["maxgrad",]
  # print(maxgrad)
   stock <- stock + fit0
   return(list(stock=stock,covars=covars))

}

#--------------------------------------------------------------------
getFmod <- function(stk, model="interaction", dfm=c(2/3, 1/2)){
  dis <- dims(stk)
  KY=floor(dfm[1] * dis$year)
  KA=ceiling(dfm[2] *dis$age)
  if(model=="interaction"){
    frm <- substitute(~ te(age, year, k=c(KA, KY)), list(KA=KA, KY=KY))
  } else if (model=="separable"){
    frm <- substitute(~ s(age, k = KA) + s(year, k = KY), list(KA=KA, KY=KY))
  }
  as.formula(frm)
}

#--------------------------------------------------------------------
getQmod <- function(idx){
  lds <- lapply(idx, dims)
  lds <- lapply(lds, function(x){
    if(x$age<=3){
      frm <- ~factor(age)
    } else {
      frm <- substitute(~s(age, k=KA), list(KA=ceiling(3/4 * x$age)))
    } 
    as.formula(frm)
  })
  lds		
}



