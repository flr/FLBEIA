###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.indices.data
# NOTE #1:      Return a list with FLIndex objects
###############################################################################
#-------------------------------------------------------------------------
#' 
#' FLBEIA easy conditioning: indices argument creator
#' 
#' create.indices.data function creates an FLIndices object
#' 
#' @param   ni Number of iterations (number).
#' @param   ns Number of seasons (number).
#' @param   yrs A vector with c(first.yr,proj.yr, last.yr) where
#'\itemize{
#'      \item first.yr: First year of simulation (number).
#'      \item proj.yr: First year of projection (number).
#'      \item last.yr: Last year of projection (number).}
#' @param   stks.data A list with the names of the stocks with indices and the following elements:
#'\itemize{
#'      \item  stk.unit: Number of units of the stock (number). 
#'      \item  stk.age.min: Minimum age class of the stock (number).
#'      \item  stk.age.max: Maximum age class of the stock (number).
#'      \item  stk_indices: Name of indices for the stock 'stk' (vector).
#'      \item  stk_ind_index.flq: Historical index data for index 'ind' of stock 'stk' (FLQuant).}
#' Optionals:
#'\itemize{
#'      \item  stk_ind_type: Type of index (character).
#'      \item  stk_ind_distribution: Name of the statistical distribution of the 'ind' index values for stock 'stk' (character).
#'      \item  stk_ind_index.var.flq: Variability in 'ind' index of stock 'stk' (FLQuant).
#'      \item  stk_ind_index.q.flq: Catchability for 'ind' index of stock 'stk' (FLQuant).
#'      \item  stk_ind_catch.n.flq: Catch at age in numbers for 'ind' index of stock 'stk'(FLQuant).
#'      \item  stk_ind_catch.wt.flq: Mean weight at age in the catch for 'ind' index of stock 'stk' (FLQuant).
#'      \item  stk_ind_effort.flq: Effort for 'ind' index of stock 'stk' (FLQuant).
#'      \item  stk_ind_sel.pattern.flq: Selection pattern for 'ind' index of stock 'stk' (FLQuant).
#'      \item  stk_ind_range.min: Minimum age in 'ind' index of stock 'stk' (number).
#'      \item  stk_ind_range.max: Maximum age in 'ind' index of stock 'stk' (number).
#'      \item  stk_ind_range.plusgroup: Plusgroup age in 'ind' index of stock 'stk' (number).
#'      \item  stk_ind_range.minyear: First year with 'ind' index data of stock 'stk' (number).
#'      \item  stk_ind_range.maxyear: Last year with 'ind' index data of stock 'stk' (number).
#'      \item  stk_ind_range.startf: Minimum age for calculating average fishing mortality for 'ind' index of stock 'stk' (number).
#'      \item  stk_ind_range.endf:   Maximum age for calculating average fishing mortality for 'ind' index of stock 'stk' (number).}
#'      
#' @return An FLIndices object.
#' 
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:        Is there any index data?
#   Section 2:        Create indices
#         2.1              Historical data per stock
#         2.2              Check dimensions
#-------------------------------------------------------------------------------

create.indices.data <- function(yrs,ns,ni,stks.data){

  #==============================================================================
  #   Section 1:        Is there any index data?
  #==============================================================================

  ind <- unlist(sapply(stks.data,function(x) grep(x, pattern="_indices",value=TRUE)))
  nms.stks.index <- unique(sub('.*?^(.*?)_indices*', '\\1', ind))
  indices <- NULL
  
  if(length(nms.stks.index)==0){
    cat("Warning:indices is NULL  \n")
  }else{

    #==============================================================================
    #   Section 2:        indices
    #==============================================================================
    
    n.stks.index   <- length(nms.stks.index)
    indices        <- vector('list',n.stks.index)
    names(indices) <- nms.stks.index
    first.yr <- yrs[["first.yr"]]
    proj.yr  <- yrs[["proj.yr"]]
    last.yr  <- yrs[["last.yr"]]
    proj.yrs       <- as.character(proj.yr:last.yr)
    hist.yrs       <- as.character(first.yr:(proj.yr-1))
    ny <- length(first.yr:last.yr)
    list.stks.unit <- lapply(stks.data, function(ch) grep(pattern="unit", ch, value = TRUE))
    list.stks.age <- lapply(stks.data, function(ch) grep(pattern="age", ch, value = TRUE))
    list.stks.flqa <- create.list.stks.flqa(nms.stks.index,yrs,ni,ns,list.stks.unit,list.stks.age)  
    list.stks.flq  <- create.list.stks.flq(nms.stks.index,yrs,ni,ns,list.stks.unit)  
    
    #-----------------------------------------------------------------------------
    #   Section 2.1:      Historical data per stock
    #-----------------------------------------------------------------------------
    
    for (i in 1:n.stks.index){
      
      stk <- nms.stks.index[i]
      cat('=============', stk,'index','=============\n')
      nms.stk.index    <- get(grep(stks.data[[stk]],pattern="_indices", value = TRUE))
      
      list.stk.index <- list()
      
      for(j in 1:length(nms.stk.index)){
        
        nm.stk.index <-  nms.stk.index[j]
        
        flqa <- list.stks.flqa[[stk]]
        flq  <- list.stks.flq[[stk]]
        
        stk.index              <- get(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_index.flq',sep=''), value = TRUE))
        stk.index.type         <- get(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_type',sep=''), value = TRUE))
        if(length(stk.index.type)==0) stk.index.type  <- NA    
        stk.index.distribution <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_distribution',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.distribution)==0) stk.index.distribution  <- NA    
        stk.index.var         <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_index.var.flq',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.var)==0) stk.index.var  <- NA    
        stk.index.q          <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_index.q.flq',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.q)==0) stk.index.q  <- NA    
        stk.index.catch.n    <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_catch.n.flq',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.catch.n)==0) stk.index.catch.n  <- NA    
        stk.index.catch.wt   <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_catch.wt.flq',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.catch.wt)==0) stk.index.catch.wt  <- NA    
        stk.index.sel.pattern<- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_sel.pattern.flq',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.sel.pattern)==0) stk.index.sel.pattern  <- NA    
        stk.index.effort     <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_effort.flq',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.effort)==0) stk.index.effort  <- NA    
        stk.index.range.min  <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.min',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.min)==0) stk.index.range.min  <- NA    
        stk.index.range.max  <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.max',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.max)==0) stk.index.range.max  <- NA    
        stk.index.range.startf    <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.startf',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.startf)==0) stk.index.range.startf  <- NA    
        stk.index.range.endf      <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.endf',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.endf)==0) stk.index.range.endf  <- NA    
        stk.index.range.plusgroup     <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.plusgroup',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.plusgroup)==0) stk.index.range.plusgroup  <- NA    
        stk.index.range.minyear     <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.minyear',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.minyear)==0) stk.index.range.minyear  <- NA    
        stk.index.range.maxyear     <- mget(grep(stks.data[[stk]],pattern=paste(stk,'_',nm.stk.index,'_range.maxyear',sep=''), value = TRUE),envir=as.environment(1))
        if(length(stk.index.range.maxyear)==0) stk.index.range.maxyear  <- NA    
        stk.unit     <- mget(grep(stks.data[[stk]],pattern=paste(stk,'.unit',sep=""), value = TRUE),envir=as.environment(1))
        if(length(stk.unit)==0) stk.unit  <- NA    
 
        #-----------------------------------------------------------------------------
        #   Section 2.2:      Check dimensions
        #-----------------------------------------------------------------------------
 
        if(stk.index.type=="FLIndex"){
          stk.flindex<- FLIndex(name = nm.stk.index, index = flqa)
        }else{
          stk.flindex<- FLIndexBiomass(name = nm.stk.index, index = flq)
          flqa <- flq}  

        if(!all(is.na(stk.index.var))){
          stk.index.var <- stk.index.var[[1]]}
        if(!all(is.na(stk.index.q))){
          stk.index.q <- stk.index.q[[1]]}
        if(!all(is.na(stk.index.effort))){
          stk.index.effort <- stk.index.effort[[1]]}
        if(!all(is.na(stk.index.catch.n))){
          stk.index.catch.n <- stk.index.catch.n[[1]]}
        if(!all(is.na(stk.index.catch.wt))){
          stk.index.catch.wt <- stk.index.catch.wt[[1]]}
        if(!all(is.na(stk.index.sel.pattern))){
          stk.index.sel.pattern <- stk.index.sel.pattern[[1]]}
        if(!all(is.na(stk.index.range.min))){
          stk.index.range.min <- stk.index.range.min[[1]]}
        if(!all(is.na(stk.index.range.max))){
          stk.index.range.max <- stk.index.range.max[[1]]}
        if(!all(is.na(stk.index.range.startf))){
          stk.index.range.startf <- stk.index.range.startf[[1]]}
        if(!all(is.na(stk.index.range.endf))){
          stk.index.range.endf <- stk.index.range.endf[[1]]}
        if(!all(is.na(stk.index.range.plusgroup))){
          stk.index.range.plusgroup <- stk.index.range.plusgroup[[1]]}
        if(!all(is.na(stk.index.range.minyear))){
          stk.index.range.minyear <- stk.index.range.minyear[[1]]}
        if(!all(is.na(stk.index.range.maxyear))){
          stk.index.range.maxyear <- stk.index.range.maxyear[[1]]}

        
        log.dim <- equal.flq.Dimnames(lflq=list(stk.index,stk.flindex@index),1:2)
        
        if(!log.dim)stop('in indices dimension \n')
        
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
      
      list.stk.index[[j]] <- stk.flindex
      }
      names(list.stk.index) <- nms.stk.index
    indices[[stk]]<- FLIndices(list.stk.index)
  }
  return(indices)
  }
}
  