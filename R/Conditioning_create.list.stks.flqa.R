###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.list.stks.flqa
# NOTE #1:      Create a list with an FLQuant for each stock with age from 
#               min.age to max.age
###############################################################################
#-------------------------------------------------------------------------
#  inputs: 
#
#   stks:     Name of all the stocks (vector)
#   first.yr: First year of simulation (number)
#   proj.yr:	First year of projection (number)
#   last.yr:	Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:	      Number of seasons (number)
#   stk.unit	Number of units of the stock (number) 
#   stk.min.age: Minimum age of the stock (number)
#   stk.max.age: Maximum age of the stock (number)
#-------------------------------------------------------------------------


create.list.stks.flqa <- function(){
  
  n.stk   <- length(stks)    
  nmy     <- as.character(first.yr:last.yr)
  list.stks.flqa <- NULL
  
  for(i in 1:n.stk){
    nmstk         <- stks[i]
    stk.age.min <- get(paste(nmstk,'.age.min',sep=''))
    stk.age.max <- get(paste(nmstk,'.age.max',sep=''))
    stk.unit    <- get(paste(nmstk,'.unit',sep=""))
    nmu <- 1:stk.unit
    nms <- 1:ns
    if(stk.unit==1) nmu <- 'unique'
    if(ns==1) nms <- 'all'
    
    if(stk.age.max-stk.age.min==0){
      stk.flqa <- FLQuant(dimnames = list(age = 'all', year = nmy, unit = nmu, 
                                          season = nms, iter = 1:ni))         
    } else{
      stk.flqa <- FLQuant(dimnames = list(age = stk.age.min:stk.age.max, year = nmy, unit = nmu, 
                                          season = nms, iter = 1:ni)) 
    }
    assign(paste(nmstk,'.flqa',sep=""),stk.flqa)
    list.stks.flqa[[i]]<- get(paste(nmstk,'.flqa',sep=""))
  }
  
  names(list.stks.flqa)<- stks
  return(list.stks.flqa) 
}
