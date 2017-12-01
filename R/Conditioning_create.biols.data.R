###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.biols.data
# NOTE #1:      Return FLBiols object called biols
###############################################################################
#-------------------------------------------------------------------------------
#
#' FLBEIA easy conditioning: biols argument creator
#' 
#' create.biols.data function creates an FLBiols object.
#'
#' @param   ni Number of iterations (number).
#' @param   ns	      Number of seasons (number).
#' @param   yrs A vector with c(first.yr,proj.yr, last.yr) where
#'\itemize{
#'      \item first.yr: First year of simulation (number).
#'      \item proj.yr:  First year of projection (number).
#'      \item last.yr:  Last year of projection (number).}
#' @param   stks.data A list with the name of the stks and the following elements:
#'\itemize{
#'      \item  stk.unit: Number of units of the stock (number). 
#'      \item  stk.age.min: Minimum age class of the stock (number).
#'      \item  stk.age.max: Maximum age class of the stock (number).
#'      \item  stk_n.flq: Numbers at age in the population(FLQuant).
#'      \item  stk_wt.flq: Weight at age of an individual (FLQuant).
#'      \item  stk_m.flq: Mortality rate at age of the population (FLQuant).	
#'      \item  stk_fec.flq: Fecundity at age	(FLQuant).	
#'      \item  stk_mat.flq: Percentage of mature individuals at age	(FLQuant).	
#'      \item  stk_spwn.flq: Proportion of time step at which spawning ocurrs	(FLQuant).
#'      \item  stk.range.plusgroup: Plusgroup age (number).
#'      \item  stk.range.minyear: Minimum year (number).
#'      \item  stk.range.maxyear: Maximum year (number).
#'      \item  stk_range.minfbar: Minimum age to calculate average fishing mortality (number).
#'      \item  stk_range.maxfbar: Maximum age to calculate average fishing mortality (number).
#'      \item  stk_biol.proj.avg.yrs:	Historic years to calculate the average of spwn, fec, m and wt for the projection (vector).}
#'      
#' @return An FLBiol object
#
#   Required functions: Create.list.stks.flqa	function
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:      Create FLBIOL per stock
#           1.1:      Historical data
#           1.2:      Projection
#   Section 2:      Create FLBiols object ('biols') with all the stocks
#   Section 3:      Return biols
#-------------------------------------------------------------------------------

create.biols.data <- function(yrs,ns,ni,stks.data){
  
  stks           <- names(stks.data)
  n.stk          <- length(stks)    
  first.yr <- yrs[["first.yr"]]
  proj.yr  <- yrs[["proj.yr"]]
  last.yr  <- yrs[["last.yr"]]
  proj.yrs       <- as.character(proj.yr:last.yr)
  hist.yrs       <- as.character(first.yr:(proj.yr-1))

  list.stks.unit <- lapply(stks.data, function(ch) grep(pattern="unit", ch, value = TRUE))
  list.stks.age <- lapply(stks.data, function(ch) grep(pattern="age", ch, value = TRUE))
  list.stks.flqa <- create.list.stks.flqa(stks,yrs,ni,ns,list.stks.unit,list.stks.age)  
  #==============================================================================
  #   Section 1:     Create FLBIOL per stock
  #==============================================================================
  list.FLBiol <- list()
  for (i in 1:n.stk){   

    nmstk        <- stks[i]
    stk.flqa <- list.stks.flqa [[nmstk]]
    
    cat('=============', nmstk,'biol','=============\n')

    #------------------------------------------------------------------------------
    #   Section 1.1:      Historical data
    #------------------------------------------------------------------------------
    stk.unit  <- get(grep(stks.data[[nmstk]],pattern=".unit", value = TRUE))
    stk.wt     <- get(grep(stks.data[[nmstk]],pattern="_wt.flq", value = TRUE))
    stk.n      <- get(grep(stks.data[[nmstk]],pattern="_n.flq", value = TRUE))
    stk.m      <- get(grep(stks.data[[nmstk]],pattern="_m.flq", value = TRUE))
    stk.fec    <- get(grep(stks.data[[nmstk]],pattern="_fec.flq", value = TRUE))
    stk.mat    <- get(grep(stks.data[[nmstk]],pattern="_mat.flq", value = TRUE))
    stk.spwn   <- get(grep(stks.data[[nmstk]],pattern="_spwn.flq", value = TRUE))
    stk.range.min       <- get(grep(stks.data[[nmstk]],pattern=".age.min", value = TRUE)) 
    stk.range.max       <- get(grep(stks.data[[nmstk]],pattern=".age.max", value = TRUE)) 
    stk.range.plusgroup <- get(grep(stks.data[[nmstk]],pattern="_range.plusgroup", value = TRUE)) 
    stk.range.minyear   <- get(grep(stks.data[[nmstk]],pattern="_range.minyear", value = TRUE)) 
    stk.range.maxyear   <- last.yr    
    stk.range.minfbar   <-  get(grep(stks.data[[nmstk]],pattern="_range.minfbar", value = TRUE)) 
    stk.range.maxfbar   <- get(grep(stks.data[[nmstk]],pattern="_range.maxfbar", value = TRUE)) 
    stk.proj.avg.yrs    <- get(grep(stks.data[[nmstk]],pattern="_biol.proj.avg.yrs", value = TRUE)) 
    stk.proj.avg.yrs    <- as.character(stk.proj.avg.yrs) 
    
    # Check the dimension names of age and years
    log.dim <- equal.flq.Dimnames(lflq=list(stk.wt,stk.n,stk.m,stk.fec,stk.spwn,stk.mat,
                                            stk.flqa[,hist.yrs,]),1:2)
    
    if(!log.dim){stop('In the dimension names of FLQuants age or years')}
    if(!(any(dim(stk.wt)[3]==c(1,stk.unit)) & any(dim(stk.n)[3]==c(1,stk.unit)) & any(dim(stk.m)[3]==c(1,stk.unit)) & 
           any(dim(stk.fec)[3]==c(1,stk.unit)) & any(dim(stk.mat)[3]==c(1,stk.unit)) & any(dim(stk.spwn)[3]==c(1,stk.unit)))){stop('Number of stock units 1 or stk.unit')}
    if(!(any(dim(stk.wt)[4]==c(1,ns)) & any(dim(stk.mat)[3]==c(1,stk.unit)) & any(dim(stk.n)[4]==c(1,ns)) & any(dim(stk.m)[4]==c(1,ns)) & 
           any(dim(stk.fec)[4]==c(1,ns)) & any(dim(stk.mat)[3]==c(1,stk.unit)) & any(dim(stk.spwn)[4]==c(1,ns)))){stop('Number of seasons 1 or ns')}
    if(!(any(dim(stk.wt)[6]==c(1,ni)) & any(dim(stk.n)[6]==c(1,ni)) & any(dim(stk.m)[6]==c(1,ni)) & 
           any(dim(stk.fec)[6]==c(1,ni)) & any(dim(stk.mat)[3]==c(1,stk.unit)) & any(dim(stk.spwn)[6]==c(1,ni)))){stop('Number of iterations 1 or ni')}
    
    # Historical NA-s transformed in 0-s
    
    stk.wt[is.na(stk.wt)]     <- 0  
    stk.n[is.na(stk.n)]       <- 0  
    stk.m[is.na(stk.m)]       <- 0  
    stk.fec[is.na(stk.fec)]   <- 0  
    stk.mat[is.na(stk.mat)]   <- 0  
    stk.spwn[is.na(stk.spwn)] <- 0  

    stk.biol   <- FLBiol(n = stk.flqa, m=stk.flqa, wt=stk.flqa, spwn= stk.flqa, name=nmstk)
    units(stk.biol) <- list(n=units(stk.n), m=units(stk.m), wt=units(stk.wt), spwn=units(stk.spwn))
    
    fec(stk.biol) <- stk.flqa
    stk.biol@n[,hist.yrs]   <- stk.n
    stk.biol@m[,hist.yrs]   <- stk.m
    stk.biol@wt[,hist.yrs]  <- stk.wt
    fec(stk.biol)[,hist.yrs] <- stk.fec
    mat(stk.biol)[,hist.yrs] <- stk.mat
    spwn(stk.biol)[,hist.yrs]<- stk.spwn
    
    stk.biol@range[1] <- stk.range.min
    stk.biol@range[2] <- stk.range.max
    stk.biol@range[3] <- stk.range.plusgroup
    stk.biol@range[4] <- stk.range.minyear
    stk.biol@range[5] <- stk.range.maxyear 
    stk.biol@range[6] <- stk.range.minfbar   
    stk.biol@range[7] <- stk.range.maxfbar  
    
    names(stk.biol@range)[6] <- 'minfbar'
    names(stk.biol@range)[7] <- 'maxfbar'
  
    #------------------------------------------------------------------------------
    #   Section 1.2:      Projection
    #------------------------------------------------------------------------------
    
    for(yr in proj.yrs){
      stk.biol@wt[,yr]  <- yearMeans(stk.wt[,stk.proj.avg.yrs])              
      fec(stk.biol)[,yr] <- yearMeans(stk.fec[,stk.proj.avg.yrs]) 
      mat(stk.biol)[,yr] <- yearMeans(stk.mat[,stk.proj.avg.yrs]) 
      stk.biol@m[,yr]   <- yearMeans(stk.m[,stk.proj.avg.yrs])            
      spwn(stk.biol)[,yr]<- yearMeans(stk.spwn[,stk.proj.avg.yrs])          
    }
    list.FLBiol[[i]] <- stk.biol
  }
  
  #==============================================================================
  #   Section 2:     Create FLBiols object ('biols') with all the stocks
  #==============================================================================
  
  names(list.FLBiol) <- stks
  biols <- FLBiols(list.FLBiol)
  
  #==============================================================================
  #   SECTION 3:     Return
  #==============================================================================
  
  return(biols)
}
