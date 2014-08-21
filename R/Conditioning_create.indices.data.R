###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.indices.data
# NOTE #1:      Return a list with FLIndex objects
###############################################################################
#-------------------------------------------------------------------------
#  inputs: 
#
#   Required:
#   first.yr:   First year of simulation (number)
#   proj.yr:    First year of projection (number)
#   last.yr:    Last year of projection (number)
#   ni:         Number of iterations (number)
#   ns:         Number of seasons (number)
#   stks_index: Names of the stocks with index (vector)
#   stk_indices:          Name of index for the stock 'stk'
#   stk_ind_index.flq:	  Historical index data of stock 'stk'(FLQuant)
#   stk_ind_index.var.flq:Variability in index of stock 'stk' (FLQuant)
#   stk_ind_index.q.flq:  Catchability by index of stock 'stk' (FLQuant)
#   stk_ind_catch.n.flq:  Number of catch at age of stock 'stk'(FLQuant)
#   stk_ind_catch.wt.flq:  Weight of catch at age of stock 'stk'(FLQuant)
#   stk_ind_effort.flq:	  Effort by index of stock 'stk'(FLQuant)
#   stk_ind_sel.pattern.flq: Selection pattern of index of stock 'stk' (FLQuant)
#   stk_ind_range.min:    Minimum age catch by index of stock 'stk'(number)
#   stk_ind_range.max:    Maximum age catch by index of stock 'stk'(number)
#   stk_ind_range.minyear: First year with index data of stock 'stk' (number)
#   stk_ind_range.maxyear: Last year with index data of stock 'stk'(number)
#   stk_ind_range.startf:	Minimum age to take into account in 'f'
#   stk_ind_range.endf:   Maximum age to take into account in 'f'
#
#   Optional:
#   stk_distribution:  Name of the stock 'stk' distribution (character)
#   stk_type:		Type of stock 'stk' (character)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        Is there any index data?
#   Section 2:        Create indices
#         2.1              Historical data per stock
#         2.2              Check dimensions
#-------------------------------------------------------------------------------

create.indices.data <- function(){

  #==============================================================================
  #   Section 1:        Is there any index data?
  #==============================================================================
  
  nms.stks.index <- 'stks_index'
  
  if(!exists(nms.stks.index)){
    indices <- NULL
    
  }else{
    
    #==============================================================================
    #   Section 2:        indices
    #==============================================================================
    
    nms.stks.index <- get(nms.stks.index)
    n.stks.index   <- length(nms.stks.index)
    indices        <- vector('list',n.stks.index)
    names(indices) <- nms.stks.index
    hist.yrs       <- as.character(first.yr:last.yr)
    
    #-----------------------------------------------------------------------------
    #   Section 2.1:      Historical data per stock
    #-----------------------------------------------------------------------------
    
    for (i in 1:n.stks.index){
      
      stk <- nms.stks.index[i]
      
      cat('=============', stk,'index','=============\n')
      
      nms.stk.index <- get(paste(stk,'_indices',sep=''))
      
      for(j in 1:length(nms.stk.index)){
        
        nm.stk.index <-  nms.stk.index[j]
        
        flqa <- create.list.stks.flqa ()[[stk]]
        flq  <- create.list.stks.flq ()[[stk]]
        
        stk.index              <- get(paste(stk,'_',nm.stk.index,'_index.flq',sep=''))                                                                
        stk.index.type         <- mget(paste(stk,'_',nm.stk.index,'_type',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]                
        stk.index.distribution <- mget(paste(stk,'_',nm.stk.index,'_distribution',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]        
        stk.index.var          <- get(paste(stk,'_',nm.stk.index,'_index.var.flq',sep=''))      
        stk.index.q            <- get(paste(stk,'_',nm.stk.index,'_index.q.flq',sep=''))        
        stk.index.catch.n      <- mget(paste(stk,'_',nm.stk.index,'_catch.n.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]   
        stk.index.catch.wt     <- mget(paste(stk,'_',nm.stk.index,'_catch.wt.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]        
        stk.index.effort       <- mget(paste(stk,'_',nm.stk.index,'_effort.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]   
        stk.index.sel.pattern  <- mget(paste(stk,'_',nm.stk.index,'_sel.pattern.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]     
      
        stk.index.range.min      <- get(paste(stk,'_',nm.stk.index,'_range.min',sep=''))    
        stk.index.range.max      <- get(paste(stk,'_',nm.stk.index,'_range.max',sep=''))    
        stk.index.range.minyear  <- mget(paste(stk,'_',nm.stk.index,'_range.minyear',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]     
        stk.index.range.maxyear  <- mget(paste(stk,'_',nm.stk.index,'_range.maxyear',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]     
        stk.index.range.startf   <- get(paste(stk,'_',nm.stk.index,'_range.startf',sep=''))
        stk.index.range.endf     <- get(paste(stk,'_',nm.stk.index,'_range.endf',sep=''))
        stk.index.range.plusgroup<- mget(paste(stk,'_',nm.stk.index,'_range.plusgroup',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]    
  
        stk.unit     <- get(paste(stk,'.unit',sep=""))
      
        #-----------------------------------------------------------------------------
        #   Section 2.2:      Check dimensions
        #-----------------------------------------------------------------------------
        
        log.dim <- equal.flq.Dimnames(lflq=list(flqa[,hist.yrs],stk.index[,hist.yrs], stk.index.var[,hist.yrs],
                                    stk.index.q[,hist.yrs]),1:2)
        if(!log.dim)stop('in indices dimension \n')
      
        if(!all(is.na(stk.index))){
          if(!(any(dim(stk.index)[3]==c(1,stk.unit))))stop('in stk.index number of units 1 or stk.unit')
          if(!(any(dim(stk.index)[4]==c(1,ns))))stop('in stk.index number of seasons 1 or ns')
          if(!(any(dim(stk.index)[6]==c(1,ni))))stop('in stk.index number of iterations 1 or ni')}
      
        if(!all(is.na(stk.index.var))){
          if(!(any(dim(stk.index.var)[3]==c(1,stk.unit))))stop('in stk.index.var number of units 1 or stk.unit')
          if(!(any(dim(stk.index.var)[4]==c(1,ns))))stop('in stk.index.var number of seasons 1 or ns')
          if(!(any(dim(stk.index.var)[6]==c(1,ni))))stop('in stk.index.var number of iterations 1 or ni')}
      
        if(!all(is.na(stk.index.q))){
          if(!(any(dim(stk.index.q)[3]==c(1,stk.unit))))stop('in stk.index.q number of units 1 or stk.unit')
          if(!(any(dim(stk.index.q)[4]==c(1,ns))))stop('in stk.index.q number of seasons 1 or ns')
          if(!(any(dim(stk.index.q)[6]==c(1,ni))))stop('in stk.index.q number of iterations 1 or ni')}

        if(!all(is.na(stk.index.effort))){
          log.dim <- equal.flq.Dimnames(lflq=list(flq[,hist.yrs],stk.index.effort[,hist.yrs]),1:2)
          if(!log.dim)stop('in stk.index.effort dimension \n')          
          if(!(any(dim(stk.index.effort)[3]==1)))stop('in stk.index.effort number of units 1 or stk.unit')
          if(!(any(dim(stk.index.effort)[4]==c(1,ns))))stop('in stk.index.effort number of seasons 1 or ns')
          if(!(any(dim(stk.index.effort)[6]==c(1,ni))))stop('in stk.index.effort number of iterations 1 or ni')}
      
        if(!all(is.na(stk.index.catch.n))){
          log.dim <- equal.flq.Dimnames(lflq=list(flqa[,hist.yrs],stk.index.catch.n[,hist.yrs]),1:2)
          if(!log.dim)stop('in stk.index.catch.n dimension \n')          
          if(!(any(dim(stk.index.catch.n)[3]==c(1,stk.unit))))stop('in stk.index.catch.n number of units 1 or stk.unit')
          if(!(any(dim(stk.index.catch.n)[4]==c(1,ns))))stop('in stk.index.catch.n number of seasons 1 or ns')
          if(!(any(dim(stk.index.catch.n)[6]==c(1,ni))))stop('in stk.index.catch.n number of iterations 1 or ni')}
      
        if(!all(is.na(stk.index.catch.wt))){
          log.dim <- equal.flq.Dimnames(lflq=list(flqa[,hist.yrs],stk.index.catch.wt[,hist.yrs]),1:2)
          if(!log.dim)stop('in stk.index.catch.wt dimension \n')          
          if(!(any(dim(stk.index.catch.wt)[3]==c(1,stk.unit))))stop('in stk.index.catch.wt number of units 1 or stk.unit')
          if(!(any(dim(stk.index.catch.wt)[4]==c(1,ns))))stop('in stk.index.catch.wt number of seasons 1 or ns')
          if(!(any(dim(stk.index.catch.wt)[6]==c(1,ni))))stop('in stk.index.catch.wt number of iterations 1 or ni')}
      
        if(!all(is.na(stk.index.sel.pattern))){
          log.dim <- equal.flq.Dimnames(lflq=list(flqa[,hist.yrs],stk.index.sel.pattern[,hist.yrs]),1:2)
          if(!log.dim)stop('in stk.index.sel.pattern dimension \n')          
          if(!(any(dim(stk.index.sel.pattern)[3]==c(1,stk.unit))))stop('in sstk.index.sel.pattern number of units 1 or stk.unit')
          if(!(any(dim(stk.index.sel.pattern)[4]==c(1,ns))))stop('in stk.index.sel.pattern number of seasons 1 or ns')
          if(!(any(dim(stk.index.sel.pattern)[6]==c(1,ni))))stop('in stk.index.sel.pattern number of iterations 1 or ni')}
        
        stk.flindex<- FLIndex(name = nm.stk.index, index = flqa)
        
        stk.flindex@catch.n[]     <- stk.index.catch.n
        stk.flindex@effort[]      <- stk.index.effort
        stk.flindex@index.q[]     <- stk.index.q
        stk.flindex@index[]       <- stk.index
        stk.flindex@index.var[]   <- stk.index.var
        stk.flindex@sel.pattern[] <- stk.index.sel.pattern 
        stk.flindex@catch.wt[]    <- stk.index.catch.wt
        stk.flindex@range[['min']]     <- stk.index.range.min
        stk.flindex@range[['max']]     <- stk.index.range.max
        stk.flindex@range[['minyear']] <- stk.index.range.minyear
        stk.flindex@range[['maxyear']] <- stk.index.range.maxyear
        stk.flindex@range[['startf']]  <- stk.index.range.startf
        stk.flindex@range[['endf']]    <- stk.index.range.endf
        stk.flindex@range[['plusgroup']]    <- stk.index.range.plusgroup
      
      assign(paste(nm.stk.index, ".FLIndex", sep = ""), stk.flindex)
    }
    indices[[stk]]<- FLIndices(sapply(paste(nms.stk.index, ".FLIndex", sep = ""), FUN = get, 
                                      envir = sys.frame(which = -1)))
    names(indices[[stk]]) <- nms.stk.index
  }
  return(indices)
  }
}
  