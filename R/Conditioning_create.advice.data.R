###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.advice.data
# NOTE #1:      Return a list with advice for each stock
###############################################################################
#-------------------------------------------------------------------------
#' 
#' FLBEIA easy conditioning: advice argument creator
#' 
#' create.advice.data function creates a list (elements: TAC, TAE and quota.share)
#' 
#' @param   ni Number of iterations (number).
#' @param   ns Number of seasons (number).
#' @param   yrs A vector with c(first.yr,proj.yr, last.yr) where
#'\itemize{
#'      \item first.yr: First year of simulation (number).
#'      \item proj.yr: First year of projection (number).
#'      \item last.yr: Last year of projection (number).}
#' @param   stks.data A list with the names of the stocks and the following elements:
#' Optionals:
#'\itemize{
#'      \item  stk_advice.TAC.flq: TAC of the stock 'stk' (FLQuant). 
#'      \item  stk_advice.TAE.flq: TAE of the stock 'stk' (FLQuant).
#'      \item  stk_advice.quota.share.flq: Quota share of the stock 'stk' (FLQuant).
#'      \item  stk_advice.avg.yrs: Historic years to calculate the average of TAC, TAE or quota share of the stock 'stk' (FLQuant).}
#'      
#' @param fleets Optional argument only required if stk_advice.quota.share is not specified. It could be the output of create_fleets_FLBEIA function (FLFleets). 
#'      
#' @return A list with TAC, TAE and quota.share elements.
#' 
#------------------------------------------------------------------------------- 

#-------------------------------------------------------------------------------
#   Section 1:        creating a list 
#   Section 2:        Historic data
#           2.1         If there is not projection data, then average
#   Section 3:        Return advice
#-------------------------------------------------------------------------------

create.advice.data<- function(yrs,ns,ni,stks.data,fleets){
  
  #==============================================================================
  #   Section 1:        creating a list 
  #==============================================================================
  
  advice        <- vector('list',3)
  names(advice) <- c('TAC','TAE','quota.share')
  
  first.yr <- yrs[["first.yr"]]
  proj.yr  <- yrs[["proj.yr"]]
  last.yr  <- yrs[["last.yr"]]
  proj.yrs       <- as.character(proj.yr:last.yr)
  hist.yrs       <- as.character(first.yr:(proj.yr-1))
  yrs <- as.character(first.yr:last.yr)
  ny <- length(first.yr:last.yr)
  
  stks <- names(stks.data)
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

  stk.advice.TAC     <- mget(grep(stks.data[[stk]],pattern="_advice.TAC.flq", value = TRUE),envir=as.environment(1))
  if(length(stk.advice.TAC)==0) stk.advice.TAC  <- NA    
  stk.advice.TAE     <- mget(grep(stks.data[[stk]],pattern="_advice.TAE.flq", value = TRUE),envir=as.environment(1))
  if(length(stk.advice.TAE)==0) stk.advice.TAE  <- NA    
  stk.advice.quota.share     <- mget(grep(stks.data[[stk]],pattern="_advice.quota.share.flq", value = TRUE),envir=as.environment(1))
  if(length(stk.advice.quota.share)==0) stk.advice.quota.share  <- NA    

  k.TAC <- 1
  if(!all(is.na(stk.advice.TAC))){
    stk.advice.TAC <- stk.advice.TAC[[1]]
    log.dim <- equal.flq.Dimnames(lflq=list(stk.advice.TAC,advice$TAC[stk,]),2)
    if(!log.dim)stop('in TAC years dimension names \n')
    if(!(any(dim(stk.advice.TAC)[4]==c(1,ns))))stop('in TAC number of seasons 1 or ns')
    if(!(any(dim(stk.advice.TAC)[6]==c(1,ni))))stop('in TAC number of iterations 1 or ni')
    if (k.TAC == 1) { units(advice$TAC) <- units(stk.advice.TAC) 
    } else { if (units(advice$TAC) != units(stk.advice.TAC)) stop('all TAC units should be equal')}
    k.TAC <- k.TAC + 1
  }
  
  k.TAE <- 1
  if(!all(is.na(stk.advice.TAE))){
    stk.advice.TAE <- stk.advice.TAE[[1]]
    log.dim <- equal.flq.Dimnames(lflq=list(stk.advice.TAE,advice$TAE[stk,]),2)
    if(!log.dim)stop('in TAE years dimension names \n')
    if(!(any(dim(stk.advice.TAE)[4]==c(1,ns))))stop('in TAE number of seasons 1 or ns')
    if(!(any(dim(stk.advice.TAE)[6]==c(1,ni))))stop('in TAE number of iterations 1 or ni')
    if (k.TAE==1) { units(advice$TAE) <- units(stk.advice.TAE) 
    } else { if (units(advice$TAE) != units(stk.advice.TAE)) stop('all TAE units should be equal')}
    k.TAE <- k.TAE + 1
  }
  
  if(!all(is.na(stk.advice.quota.share))){
    stk.advice.quota.share <- stk.advice.quota.share[[1]]
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
          if(sum(flinfo[stk,which(sapply(strsplit(colnames(flinfo), "&&"), function(x) x[1]) == f)]) == 0)  next
            # Historic period
            advice$quota.share[[stk]][f,hist.yrs,] <-  apply(catchWStock.f(fleets[[f]], stk), c(2,6), sum)[,hist.yrs,]/totC
            # Projection period
            advice$quota.share[[stk]][f,proj.yrs,] <- yearMeans(advice$quota.share[[stk]][f,stk.proj.avg.yrs,])
      }
    }
    advice$quota.share[[stk]][is.na(advice$quota.share[[stk]])] <- 0
  }
  #==============================================================================
  #   Section 3:        Return advice
  #==============================================================================
  
  return(advice)
}

