###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.biols.data
# NOTE #1:      Return FLBiols object called biols
###############################################################################
#-------------------------------------------------------------------------
#  inputs: 
#
#   first.yr: First year of simulation (number)
#   proj.yr:  First year of projection (number)
#   last.yr:	Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:	      Number of seasons (number)
#   stks:     Name of all the stocks (vector)
#   stk.unit	Number of units of the stock (number) 
#   stk.min.age: Minimum age of the stock (number)
#   stk.max.age: Maximum age of the stock (number)
#   stk_n.flq:   Abundance at age age (FLQuant)
#   stk_wt.flq:  Weight age age (FLQuant)
#   stk_m.flq:   Natural mortality age age (FLQuant)	
#   stk_fec.flq: Fecundity at age	(FLQuant)	
#   stk_spwn.flq:   Spawning at age	(FLQuant)
#   stk_range.min:	Minimum age (number)
#   stk_range.max:	Maximum age (number)
#   stk.range.minyear:  Minimum year (number)
#   stk_range.minfbar:	Minimum year to take into account in the calculation of 'f' (number)
#   stk_range.maxfbar:	Maximum year to take into account in the calculation of 'f' (number)
#   stk_biol.proj.avg.yrs:	historic years to calculate the average of spwn, fec, m and wt (vector)
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

create.biols.data <- function(){
  
  n.stk          <- length(stks)    
  proj.yrs       <- as.character(proj.yr:last.yr)
  hist.yrs       <- as.character(first.yr:(proj.yr-1))
  list.stks.flqa <- create.list.stks.flqa()  
  #==============================================================================
  #   Section 1:     Create FLBIOL per stock
  #==============================================================================
  
  for (i in 1:n.stk){   

    nmstk        <- stks[i]
    stk.flqa <- list.stks.flqa [[nmstk]]
    
    cat('=============', nmstk,'biol','=============\n')

    #------------------------------------------------------------------------------
    #   Section 1.1:      Historical data
    #------------------------------------------------------------------------------
    stk.unit  <- get(paste(nmstk,'.unit',sep=""))
    stk.wt     <- get(paste(nmstk,'_wt.flq',sep=""))
    stk.n      <- get(paste(nmstk,'_n.flq',sep=""))
    stk.m      <- get(paste(nmstk,'_m.flq',sep=""))
    stk.fec    <- get(paste(nmstk,'_fec.flq',sep=""))
    stk.spwn   <- get(paste(nmstk,'_spwn.flq',sep=""))
    stk.range.min       <- get(paste(nmstk,'_range.min',sep=""))
    stk.range.max       <- get(paste(nmstk,'_range.max',sep=""))
    stk.range.plusgroup <- get(paste(nmstk,'_range.plusgroup',sep=""))
    stk.range.minyear   <- get(paste(nmstk,'_range.minyear',sep=""))
    stk.range.maxyear   <- last.yr    
    stk.range.minfbar   <- get(paste(nmstk,'_range.minfbar',sep=""))
    stk.range.maxfbar   <- get(paste(nmstk,'_range.maxfbar',sep=""))
    stk.proj.avg.yrs    <- get(paste(nmstk,'_biol.proj.avg.yrs',sep=""))
    stk.proj.avg.yrs    <- as.character(stk.proj.avg.yrs) 
    
    # Check the dimension names of age and years
    log.dim <- equal.flq.Dimnames(lflq=list(stk.wt,stk.n,stk.m,stk.fec,stk.spwn,
                                            stk.flqa[,hist.yrs,]),1:2)
    
    if(!log.dim){stop('In the dimension names of FLQuants age or years')}
    if(!(any(dim(stk.wt)[3]==c(1,stk.unit)) & any(dim(stk.n)[3]==c(1,stk.unit)) & any(dim(stk.m)[3]==c(1,stk.unit)) & 
           any(dim(stk.fec)[3]==c(1,stk.unit)) & any(dim(stk.spwn)[3]==c(1,stk.unit)))){stop('Number of stock units 1 or stk.unit')}
    if(!(any(dim(stk.wt)[4]==c(1,ns)) & any(dim(stk.n)[4]==c(1,ns)) & any(dim(stk.m)[4]==c(1,ns)) & 
           any(dim(stk.fec)[4]==c(1,ns)) & any(dim(stk.spwn)[4]==c(1,ns)))){stop('Number of seasons 1 or ns')}
    if(!(any(dim(stk.wt)[6]==c(1,ni)) & any(dim(stk.n)[6]==c(1,ni)) & any(dim(stk.m)[6]==c(1,ni)) & 
           any(dim(stk.fec)[6]==c(1,ni)) & any(dim(stk.spwn)[6]==c(1,ni)))){stop('Number of iterations 1 or ni')}
    
    # Historical NA-s transformed in 0-s
    
    stk.wt[is.na(stk.wt)]     <- 0
    stk.n[is.na(stk.n)]       <- 0
    stk.m[is.na(stk.m)]       <- 0
    stk.fec[is.na(stk.fec)]   <- 0
    stk.spwn[is.na(stk.spwn)] <- 0
    
    stk.biol   <- FLBiol(n = stk.flqa, m=stk.flqa, wt=stk.flqa, fec=stk.flqa, spwn= stk.flqa, name=nmstk)

    stk.biol@n[,hist.yrs]   <- stk.n
    stk.biol@m[,hist.yrs]   <- stk.m
    stk.biol@wt[,hist.yrs]  <- stk.wt
    stk.biol@fec[,hist.yrs] <- stk.fec
    stk.biol@spwn[,hist.yrs]<- stk.spwn
    
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
      stk.biol@fec[,yr] <- yearMeans(stk.fec[,stk.proj.avg.yrs])          
      stk.biol@m[,yr]   <- yearMeans(stk.m[,stk.proj.avg.yrs])            
      stk.biol@spwn[,yr]<- yearMeans(stk.spwn[,stk.proj.avg.yrs])          
    }

    assign(paste(nmstk,".biol",sep=""),stk.biol)
  }
  
  #==============================================================================
  #   Section 2:     Create FLBiols object ('biols') with all the stocks
  #==============================================================================
  
  biols        <- FLBiols(sapply(paste(stks,".biol",sep=""),FUN=get, envir=sys.frame(which=-1)))
  names(biols) <- stks
    
  #==============================================================================
  #   SECTION 3:     Return
  #==============================================================================
  
  return(biols)
}
