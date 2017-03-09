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
#   yrs:      a vector with c(first.yr,proj.yr,last.yr)
#                 first.yr:	First year of simulation (number)
#                 proj.yr:	First year of projection (number)
#                 last.yr:	Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:	      Number of seasons (number)
#   list.stks.unit: a list with the name of the stks and each 
#                     stk includes the number of units 
#-------------------------------------------------------------------------


create.list.stks.flq <- function(stks,yrs,ni,ns,list.stks.unit){
  
  n.stk   <- length(stks)  
  first.yr <- yrs[["first.yr"]]
  proj.yr  <- yrs[["proj.yr"]]
  last.yr  <- yrs[["last.yr"]]
  nmy     <- as.character(first.yr:last.yr)
  list.stks.flq <- NULL
  
  for(i in 1:n.stk){
    nmstk     <- stks[i]  
    stk.unit <- get(grep(list.stks.unit[[nmstk]],pattern="unit", value = TRUE))
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
