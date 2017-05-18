###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.SRs.data
# NOTE #1:      Return a list with FLSRsim objects
###############################################################################
#-------------------------------------------------------------------------
#' 
#' FLBEIA easy conditioning: SRs argument creator
#' 
#' create.BDs.data function creates a list of FLBDsim objects
#' 
#' @param   ni Number of iterations (number).
#' @param   ns Number of seasons (number).
#' @param   yrs A vector with c(first.yr,proj.yr, last.yr) where:
#'  \itemize{
#'      \item first.yr: First year of simulation (number).
#'      \item proj.yr: First year of projection (number).
#'      \item last.yr: Last year of projection (number).}
#' @param   stks.data A list with the name of the stks and the following elements:
#'  \itemize{
#'      \item stk.unit: Number of units of the stock (number).
#'      \item  stk.age.min: Minimum age class of the stock (number).
#'      \item  stk.age.max: Maximum age class of the stock (number).
#'      \item stk_sr.model: Name of the model to simulate recruitment (character).
#'      \item stk_params.n: Number of parameters (number).
#'      \item stk_params.name: Name of the parameters (vector).
#'      \item stk_params.array:	Parameter values (array).
#'      \item stk_rec.flq: Recruitment values (FLQuant).
#'      \item stk_ssb.flq: Spawning stock biomass values (FLQuant).
#'      \item stk_proportion.flq: Recruitment distribution in each time step as a proportion (FLQuant, values between 0 and 1).
#'      \item stk_prop.avg.yrs: Historical years to calculate the proportion average (vector).
#'      \item stk_timelag.matrix: Timelag between the spawning an recruitment (matrix [2, number of seasons]). For details see FLSRsim.
#'      \item stk_range.plusgroup: Plusgroup age (number).
#'      \item stk_range.minyear: Minimum year (number).}
#'  Optionals:
#'  \itemize{
#'      \item stk_uncertainty.flq: Uncertainty (FLQuant).} 
#'   
#' @return A list of FLSRsim objects.
#'   
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        FLSRSim for each stock
#         1.1         Historic and projection data(params,uncertainty*)
#         1.2              Check dimensions
#         1.3              Projection
#   Section 2:        SRs 
#-------------------------------------------------------------------------------

create.SRs.data <- function(yrs,ns,ni,stks.data){  
  
  ind <- unlist(sapply(stks.data,function(x) grep(x, pattern="sr.model",value=TRUE)))
  nmstks <- unique(sub('.*?^(.*?)_sr.model*', '\\1', ind))
  SRs <- NULL
  
  if(length(nmstks)!=0){

    n.stk.SR     <- length(nmstks)

    first.yr <- yrs[["first.yr"]]
    proj.yr  <- yrs[["proj.yr"]]
    last.yr  <- yrs[["last.yr"]]
    proj.yrs       <- as.character(proj.yr:last.yr)
    hist.yrs       <- as.character(first.yr:(proj.yr-1))
    ny <- length(first.yr:last.yr)
    list.stks.unit <- lapply(stks.data, function(ch) grep(pattern="unit", ch, value = TRUE))
    list.stks.age <- lapply(stks.data, function(ch) grep(pattern="age", ch, value = TRUE))
    list.stks.flqa <- create.list.stks.flqa(nmstks,yrs,ni,ns,list.stks.unit,list.stks.age)  
    list.stks.flq  <- create.list.stks.flq(nmstks,yrs,ni,ns,list.stks.unit)  
        
    #==============================================================================
    #   Section 1:         Create FLSRSim for each stock
    #==============================================================================
    
    list.SRs <- list()
    for(i in 1:n.stk.SR){
      
      nmstk       <- nmstks[i]

      
      cat('=============', nmstk,'SR','=============\n')
      
      flqa.stk    <-list.stks.flqa[[nmstk]][1,,1]  #age=1 and unit=1
      flq.stk     <- list.stks.flq[[nmstk]][,,1]
      
      #-----------------------------------------------------------------------------
      #    1.1         Historic and projection data(params,uncertainty*)
      #-----------------------------------------------------------------------------

        stk.model           <- get(grep(stks.data[[nmstk]],pattern="_sr.model", value = TRUE))
        stk.unit            <- get(grep(stks.data[[nmstk]],pattern=".unit", value = TRUE))
        stk.rec             <- get(grep(stks.data[[nmstk]],pattern="_rec.flq", value = TRUE))
        stk.ssb             <- get(grep(stks.data[[nmstk]],pattern="_ssb.flq", value = TRUE))
        stk.params.n        <- get(grep(stks.data[[nmstk]],pattern="_params.n", value = TRUE))
        stk.params.name     <- get(grep(stks.data[[nmstk]],pattern="_params.name", value = TRUE))
        stk.params          <- get(grep(stks.data[[nmstk]],pattern="_params.array", value = TRUE))
        stk.proportion      <- get(grep(stks.data[[nmstk]],pattern="_proportion.flq", value = TRUE))
        stk.timelag         <- get(grep(stks.data[[nmstk]],pattern="_timelag.matrix", value = TRUE))
        stk.range.min       <- get(grep(stks.data[[nmstk]],pattern=".age.min", value = TRUE))
        stk.range.max       <- get(grep(stks.data[[nmstk]],pattern=".age.max", value = TRUE))
        stk.range.plusgroup <- get(grep(stks.data[[nmstk]],pattern="_range.plusgroup", value = TRUE))
        stk.range.minyear   <- get(grep(stks.data[[nmstk]],pattern="_range.minyear", value = TRUE))
        stk.uncertainty     <- mget(grep(stks.data[[nmstk]],pattern="_uncertainty.flq", value = TRUE),envir=as.environment(1))
        if(length(stk.uncertainty)==0) stk.uncertainty  <- NA    
       
        params <- array(dim = c(stk.params.n,ny,ns,ni), 
                 dimnames = list(param=ac(1:stk.params.n),year = ac(first.yr:last.yr),season=ac(1:ns), iter = 1:ni))
                      
        stk.sr <- FLSRsim(name = nmstk, model =stk.model, rec = flqa.stk, 
                          ssb = flq.stk,params = params, uncertainty = flqa.stk,
                          proportion = flqa.stk, covar=FLQuants())
        
        dimnames(stk.sr@params)$param <-stk.params.name
      
        stk.sr@timelag[]      <- stk.timelag
        stk.sr@range[['min']] <- stk.range.min
        stk.sr@range[['max']] <- stk.range.max
        stk.sr@range[['plusgroup']] <- stk.range.plusgroup
        stk.sr@range[['minyear']] <- stk.range.minyear
        stk.sr@range[['maxyear']] <- last.yr

      #-----------------------------------------------------------------------------
      #       1.2         Check dimensions
      #-----------------------------------------------------------------------------
      
        if(!all(is.na(stk.rec))){
          log.dim <- equal.flq.Dimnames(lflq=list(stk.rec,stk.sr@rec[,hist.yrs]),2)
          if(!log.dim)stop('in SR recruitment year dimension names \n')
          if(!(any(dim(stk.rec)[3]==c(1,stk.unit))))stop('in rec number of stock units 1 or stk.unit')
          if(!(any(dim(stk.rec)[4]==c(1,ns))))stop('in rec number of seasons 1 or ns')
          if(!(any(dim(stk.rec)[6]==c(1,ni))))stop('in rec number of iterations 1 or ni')
        }else{cat('SR recruitment all NA-s \n')}
      
        if(!all(is.na(stk.ssb))){
          log.dim <- equal.flq.Dimnames(lflq=list(stk.ssb,stk.sr@ssb[,hist.yrs]),2)
          if(!log.dim)stop('in SR ssb year dimension names \n')
          if(!(any(dim(stk.ssb)[3]==c(1,stk.unit))))stop('in ssb number of stock units 1 or stk.unit')
          if(!(any(dim(stk.ssb)[4]==c(1,ns))))stop('in ssb number of seasons 1 or ns')
          if(!(any(dim(stk.ssb)[6]==c(1,ni))))stop('in ssb number of iterations 1 or ni')
        }else{cat('SR ssb all NA-s \n')}
      
        if(!all(is.na(stk.uncertainty))){
          stk.uncertainty <- stk.uncertainty[[1]]
          log.dim <- equal.flq.Dimnames(lflq=list(stk.uncertainty,stk.sr@uncertainty),2)
          if(!log.dim)stop('in SR uncertainty year dimension names \n')
          if(!(any(dim(stk.uncertainty)[3]==c(1,stk.unit))))stop('in uncertainty number of stock units 1 or stk.unit')
          if(!(any(dim(stk.uncertainty)[4]==c(1,ns))))stop('in uncertainty number of seasons 1 or ns')
          if(!(any(dim(stk.uncertainty)[6]==c(1,ni))))stop('in uncertainty number of iterations 1 or ni')
        }else{stk.uncertainty=1
              cat('SR uncertainty = 1 \n')}
      
      if(!all(is.na(stk.proportion))){
        log.dim <- equal.flq.Dimnames(lflq=list(stk.proportion,stk.sr@proportion),2)
        if(!log.dim)stop('in SR proportion year dimension names \n')
        if(!(any(dim(stk.proportion)[3]==c(1,stk.unit))))stop('in proportion number of stock units 1 or stk.unit')
        if(!(any(dim(stk.proportion)[4]==c(1,ns))))stop('in proportion number of seasons 1 or ns')
        if(!(any(dim(stk.proportion)[6]==c(1,ni))))stop('in proportion number of iterations 1 or ni')}

      log.dim <- equal.flq.Dimnames(lflq=list(stk.params,stk.sr@params),1:3)
      if(!log.dim)stop('in SR parameter dimension names \n')
    
      stk.sr@rec[,hist.yrs] <- stk.rec
      stk.sr@ssb[,hist.yrs] <- stk.ssb
      stk.sr@params[] <- stk.params
      stk.sr@uncertainty[] <- stk.uncertainty
      stk.sr@proportion[]  <- stk.proportion
      
      #-----------------------------------------------------------------------------
      #       1.3         Projection 
      #-----------------------------------------------------------------------------
      
        if(!any(is.na(stk.params[,proj.yrs,,]))){
          if(!all(dim(stk.params)==dim(stk.sr@params))){
           stop('in SR parameters dimension names \n')}
        }else{stop('SR parameters all NA-s \n')}
    
        if(all(is.na(stk.proportion[,proj.yrs]))){
          avg.yrs   <- get(paste(nmstk,'_prop.avg.yrs',sep=""))
          for(ss in 1:ns){
          stk.sr@proportion[,proj.yrs,,ss]  <- yearMeans(stk.proportion[,avg.yrs,,ss])}
        }
        if(any(is.na(stk.sr@uncertainty[,proj.yrs]))){
          stop('Na values in uncertainty in the projection years')}

        list.SRs[[i]] <- stk.sr 
     }
    
    #==============================================================================
    #   Section 2:        Save SRs 
    #==============================================================================
    
    names(list.SRs) <- nmstks
    SRs <- list.SRs
    }  
    
  return(SRs)
  }  

