###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.SRs.data
# NOTE #1:      Return a list with FLSRsim objects
###############################################################################
#-------------------------------------------------------------------------
#  inputs: 
#
#   Required:
#   first.yr: First year of simulation (number)
#   proj.yr: First year of projection (number)
#   last.yr: Last year of projection (number)
#   ni: Number of iterations (number)
#   ns:	Number of seasons (number)
#   stk_sr.model: Name of the SR model (character)
#   stk_params.n:	Number of the parameters in the model (number)
#   stk_params.name:	Name of the parameters (vector)
#   stk_params.array:	Parameter values(array)
#   stk_rec.flq:	Recruitment values (FLQuant)
#   stk_ssb.flq:	SSB values (FLQuant)
#   stk_proportion.flq:	Proportion of recruits per season (FLQuant) 
#   stk_prop.avg.yrs:	Historical years to calculate the average of proportion (vector) 
#   stk_timelag.matrix: timelag of spawning in years and season (matrix)
#
# (optional)
#   stk_uncertainty.flq: Uncertainty (FLQuant)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        FLSRSim for each stock
#         1.1         Historic and projection data(params,uncertainty*)
#         1.2              Check dimensions
#         1.3              Projection
#   Section 2:        SRs 
#-------------------------------------------------------------------------------

create.SRs.data <- function(path){  
  
  nmstks       <- unique(sub('.*?^(.*?)_sr.model*', '\\1', ls(pattern='_sr.model',envir=sys.frame(which=0))))
  n.stk.SR     <- length(nmstks)
  
  if(n.stk.SR != 0){
    
    list.stks.flqa <- create.list.stks.flqa ()
    list.stks.flq  <- create.list.stks.flq ()
    
    ny <- length(first.yr:last.yr)
    hist.yrs <- as.character(first.yr:(proj.yr-1))
    proj.yrs <- as.character(proj.yr:last.yr)

    #==============================================================================
    #   Section 1:         Create FLSRSim for each stock
    #==============================================================================
    
    for(i in 1:n.stk.SR){
      
      nmstk       <- nmstks[i]
      
      cat('=============', nmstk,'SR','=============\n')
      
      flqa.stk    <-list.stks.flqa[[nmstk]][1,,1]  #age=1 and unit=1
   
      flq.stk     <- list.stks.flq[[nmstk]][,,1]
      
      #-----------------------------------------------------------------------------
      #    1.1         Historic and projection data(params,uncertainty*)
      #-----------------------------------------------------------------------------

        stk.model        <- get(paste(nmstk,'_sr.model',sep=''))
        stk.unit         <- get(paste(nmstk,'.unit',sep=""))
        stk.rec          <- get(paste(nmstk,'_rec.flq',sep=""))
        stk.ssb          <- get(paste(nmstk,'_ssb.flq',sep=""))
        stk.params.n     <- get(paste(nmstk,'_params.n',sep=''))
        stk.params.name  <- get(paste(nmstk,'_params.name',sep=''))
        stk.params          <- get(paste(nmstk,'_params.array',sep=""))
        stk.proportion      <- get(paste(nmstk,'_proportion.flq',sep=""))
        stk.timelag         <- get(paste(nmstk,'_timelag.matrix',sep=""))
        stk.range.min       <- get(paste(nmstk,'_range.min',sep=""))
        stk.range.max       <- get(paste(nmstk,'_range.max',sep=""))
        stk.range.plusgroup <- get(paste(nmstk,'_range.plusgroup',sep=""))
        stk.range.minyear   <- get(paste(nmstk,'_range.minyear',sep=""))
      
        stk.uncertainty    <- mget(paste(nmstk,'_uncertainty.flq',sep=''),
                                   envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
      
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
        stk.sr@range[['maxyear']] <- proj.yr-1

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
      
        assign(paste(nmstk,'.sr',sep=''),stk.sr)
    
     }
    
    #==============================================================================
    #   Section 2:        Save SRs 
    #==============================================================================
    
      all.stk.sr <- paste(nmstks,'.sr',sep='')
      SRs <-  sapply(all.stk.sr, get, envir=sys.frame(sys.parent(-1)), simplify=FALSE)
      names(SRs)<- nmstks
    }  
    
  return(SRs)
  
}
