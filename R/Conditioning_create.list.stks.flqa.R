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
#   yrs:      a vector with c(first.yr,proj.yr,last.yr)
#                 first.yr:	First year of simulation (number)
#                 proj.yr:	First year of projection (number)
#                 last.yr:	Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:	      Number of seasons (number)
#   list.stks.unit: a list with the name of the stks and each 
#                     stk includes the number of units
#   list.stks.age: a list with the name of the stks and each 
#                     stk includes min.age and max.age
#-------------------------------------------------------------------------


create.list.stks.flqa <- function(stks,yrs,ni,ns,list.stks.unit,list.stks.age){
  
  n.stk   <- length(stks)    
  first.yr <- yrs[["first.yr"]]
  proj.yr  <- yrs[["proj.yr"]]
  last.yr  <- yrs[["last.yr"]]
  nmy     <- as.character(first.yr:last.yr)
  list.stks.flqa <- NULL
  
  for(i in 1:n.stk){
    nmstk         <- stks[i]
    stk.age.min <- get(grep(list.stks.age[[nmstk]],pattern="age.min", value = TRUE))
    stk.age.max <- get(grep(list.stks.age[[nmstk]],pattern="age.max", value = TRUE))
    stk.unit <- get(grep(list.stks.unit[[nmstk]],pattern="unit", value = TRUE))
    nmu <- 1:stk.unit
    nms <- 1:ns
    if(stk.unit==1) nmu <- 'unique'
    if(ns==1) nms <- 'all'
    
    if(all(is.na(c(stk.age.max,stk.age.min)))){
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
