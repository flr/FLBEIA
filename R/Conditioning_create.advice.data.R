###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.advice.data
# NOTE #1:      Return a list with advice for each stock
###############################################################################

#-------------------------------------------------------------------------
#  inputs: 
#
#   Required:
#   first.yr: First year of simulation (number)
#   proj.yr:  First year of projection (number)
#   last.yr:  Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:       Number of seasons (number)
#
#   Optional:
#   fleets: an object called fleets, an FLFleets object, it could be the output of create_fleets_FLBEIA function. (FLFleets) 
#          Fleets is required when stk_advice.quota.share is not specify.
#   stk_advice.TAC.flq: TAC of the stock 'stk' (FLQuant)
#   stk_advice.TAE.flq:	TAE of the stock 'stk' (FLQuant)
#   stk_advice.quota.share.flq:	Quota share of the stock (FLQuant)
#   stk_advice.avg.yrs:         Historic years to calculate the average of TAC, TAE or quota share of the stock 'stk'
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        creating a list 
#   Section 2:        Historic data
#           2.1         If there is not projection data, then average
#   Section 3:        Return advice
#-------------------------------------------------------------------------------

create.advice.data<- function(){
  
  #==============================================================================
  #   Section 1:        creating a list 
  #==============================================================================
  
  advice        <- vector('list',3)
  names(advice) <- c('TAC','TAE','quota.share')
  
  yrs <- as.character(first.yr:last.yr)
  hist.yrs <- as.character(first.yr:(proj.yr-1))
  proj.yrs <- as.character(proj.yr:last.yr)
  n.stk <- length(stks)
  
  #==============================================================================
  #   Section 2:        Historic data
  #==============================================================================

  advice$TAC        <- FLQuant(dimnames = list(stock = stks, year = yrs), iter = ni)
  advice$TAE        <- FLQuant(dimnames = list(stock = stks, year = yrs), iter = ni)
  advice$quota.share <- vector('list', n.stk)

  names(advice$quota.share) <- stks
  flinfo <- stock.fleetInfo(fleets)
  flnms  <- names(fleets)
  
  for(i in 1:n.stk){
  
  stk <- stks[i]
  
  advice$quota.share[[stk]]   <- FLQuant(dimnames=list(fleet = flnms, year = yrs, iter = 1:ni))

  stk.advice.TAC    <- mget(paste(stk,'_advice.TAC.flq',sep='') ,envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
  stk.advice.TAE    <- mget(paste(stk,'_advice.TAE.flq',sep='') ,envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
  stk.advice.quota.share    <- mget(paste(stk,'_advice.quota.share.flq',sep='') ,envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]

  if(!all(is.na(stk.advice.TAC))){
    log.dim <- equal.flq.Dimnames(lflq=list(stk.advice.TAC,advice$TAC[stk,]),2)
    if(!log.dim)stop('in TAC years dimension names \n')
    if(!(any(dim(stk.advice.TAC)[4]==c(1,ns))))stop('in TAC number of seasons 1 or ns')
    if(!(any(dim(stk.advice.TAC)[6]==c(1,ni))))stop('in TAC number of iterations 1 or ni')}
  
  if(!all(is.na(stk.advice.TAE))){
    log.dim <- equal.flq.Dimnames(lflq=list(stk.advice.TAE,advice$TAE[stk,]),2)
    if(!log.dim)stop('in TAE years dimension names \n')
    if(!(any(dim(stk.advice.TAE)[4]==c(1,ns))))stop('in TAE number of seasons 1 or ns')
    if(!(any(dim(stk.advice.TAE)[6]==c(1,ni))))stop('in TAE number of iterations 1 or ni')}
  
  if(!all(is.na(stk.advice.quota.share))){
    log.dim <- equal.flq.Dimnames(lflq=list(stk.advice.quota.share,advice$quota.share[[stk]]),2)
    if(!log.dim)stop('in quota share years dimension names \n')
    if(!(any(dim(stk.advice.quota.share)[4]==c(1,ns))))stop('in quota share number of seasons 1 or ns')
    if(!(any(dim(stk.advice.quota.share)[6]==c(1,ni))))stop('in quota share number of iterations 1 or ni')}
  
  advice$TAC[stk,] <- stk.advice.TAC
  advice$TAE[stk,] <- stk.advice.TAE
  
  advice$quota.share[[stk]][] <- stk.advice.quota.share
  
  #-----------------------------------------------------------------------------
  #   2.1       If there is not projection data, then average
  #-----------------------------------------------------------------------------
 
  stk.proj.avg.yrs <- as.character(get(paste(stk,'_advice.avg.yrs',sep='')))
  
  if(any(is.na(advice$TAC[stk,proj.yrs]))){
    advice$TAC[stk,proj.yrs] <- yearMeans(advice$TAC[stk,stk.proj.avg.yrs,])
  }
  
  if(any(is.na(advice$TAE[stk,proj.yrs]))){
    advice$TAE[stk,proj.yrs] <- yearMeans(advice$TAE[stk,stk.proj.avg.yrs,])
  }
  
  if(all(is.na(advice$quota.share[[stk]][,proj.yrs]))){
    totC <- apply(catchWStock(fleets, stk), c(2,6), sum)[,hist.yrs,]
    stk.proj.avg.yrs <- as.character(get(paste(stk,'_advice.avg.yrs',sep='')))
      for(f in flnms){
          if(sum(flinfo[stk,which(sapply(strsplit(colnames(flinfo), "&&"), function(x) x[1]) == f)]) == 0)   next
            # Historic period
            advice$quota.share[[stk]][f,hist.yrs,] <-  apply(catchWStock.f(fleets[[f]], stk), c(2,6), sum)[,hist.yrs,]/totC
            # Projection period
            advice$quota.share[[stk]][f,proj.yrs,] <- yearMeans(advice$quota.share[[stk]][f,stk.proj.avg.yrs,])
      }
    }
  }
  #==============================================================================
  #   Section 3:        Return advice
  #==============================================================================
  
  return(advice)
}

