###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Create.list.stks.flq 
# NOTE #1:      Create a list with an FLQuant for each stock with age <- 'all'
###############################################################################
#-------------------------------------------------------------------------
#  inputs: 
#
#   stks:     Name of all the stocks (vector)
#   first.yr:	First year of simulation (number)
#   proj.yr:	First year of projection (number)
#   last.yr:	Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:	      Number of seasons (number)
#   stk.unit	Number of units of the stock (number)
#-------------------------------------------------------------------------


create.list.stks.flq <- function(){
  
  n.stk   <- length(stks)    
  nmy     <- as.character(first.yr:last.yr)
  list.stks.flq <- NULL
  
  for(i in 1:n.stk){
    nmstk     <- stks[i]    
    stk.unit  <- get(paste(nmstk,'.unit',sep=""))
    nmu <- 1:stk.unit
    nms <- 1:ns
    if(stk.unit==1) nmu <- 'unique'
    if(ns==1) nms <- 'all'
    
    stk.flq <- FLQuant(dimnames = list(age = 'all', year = nmy, unit = nmu, 
                                       season = nms, iter = 1:ni)) 
    
    assign(paste(nmstk,'.flq',sep=""),stk.flq)
    list.stks.flq[[i]]<- get(paste(nmstk,'.flq',sep=""))
  }

  names(list.stks.flq)<- stks
  return(list.stks.flq)  
}
