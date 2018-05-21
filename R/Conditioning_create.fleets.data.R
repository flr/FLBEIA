###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI                      
# TITLE:        create.fleets.data
# NOTE #1:      Return FLFleets object called fleets
###############################################################################
#-------------------------------------------------------------------------
#' 
#' FLBEIA easy conditioning: fleets argument creator
#' 
#' create.fleets.data function creates an FLFleetsExt object
#' 
#' @param   ni Number of iterations (number).
#' @param   ns Number of seasons (number).
#' @param   yrs A vector with c(first.yr,proj.yr, last.yr) where
#'\itemize{
#'      \item first.yr: First year of simulation (number).
#'      \item proj.yr: First year of projection (number).
#'      \item last.yr: Last year of projection (number).}
#' @param  fls.data A list with the name of the fleets and the following elements:
#'\itemize{
#'      \item fl.met: Name of the metiers in the fleet 'fl' (vector).
#'      \item fl.met.stks: Name of the stocks in the metier 'met' and fleet 'fl' (vector).
#'      \item fl_effort.flq: 'fl' fleet's effort (FLQuant).
#'      \item fl.met_effshare.flq: 'fl' fleet and 'met' metier's effort share (FLQuant).
#'      \item fl.met.stk_landings.n.flq: 'fl' fleet,'met' metier and 'stk' stock's landings in numbers at age.
#'      \item fl_proj.avg.yrs: Historic years to calculate the average of effort, fcost, crewshare, capacity in 'fl' fleet (vector) in the projection years.}
#' Optionals:
#'\itemize{
#'      \item fl_capacity.flq: 'fl' fleet's capacity (FLQuant).
#'      \item fl_fcost.flq: 'fl' fleet's fixed cost (FLQuant).
#'      \item fl_crewshare.flq: 'fl' fleet's crewshare (FLQuant).
#'      \item fl.met_vcost.flq: 'fl' fleet and 'met' metier's variable costs (FLQuant).
#'      \item fl.met.stk_landings.wt.flq: 'fl' fleet,'met' metier and 'stk' stock's mean weight of landings at age.
#'      \item fl.met.stk_discards.n.flq: 'fl' fleet,'met' metier and 'stk' stock's discards in numbers at age (FLQuant).
#'      \item fl.met.stk_discards.wt.flq: 'fl' fleet,'met' metier and 'stk' stock's mean weight of discards at age.
#'      \item fl.met.stk_price.flq: 'fl' fleet,'met' metier and 'stk' stock's price at age (FLQuant).
#'      \item fl.met.stk_alpha.flq:	'fl' fleet,'met' metier and 'stk' stock's Cobb Douglass alpha parameter (FLQuant).
#'      \item fl.met.stk_beta.flq: 'fl' fleet,'met' metier and 'stk' stock's Cobb Douglass beta parameter (FLQuant).
#'      \item fl.met.stk_catch.q.flq: 'fl' fleet,'met' metier and 'stk' stock's Cobb Douglass catch.q parameter (FLQuant).
#'      \item fl.met_proj.avg.yrs: Historic years to calculate the average of effshare,vcost for 'fl' fleet and 'met' metier in the projection period (vector).
#'      \item fl.met.stk_proj.avg.yrs: Historic years to calculate the average of landings.wt, discards.wt, landings.sel, discards.sel 
#'                                   alpha,beta,catch.q for 'fl' fleet, 'met' metier and 'stk' stock in the projection years (vector).}
#'                                   
#' @param   stks.data A list with the name of the stocks and with the next elements:
#'\itemize{
#'      \item  stk.unit: Number of units of the stock (number). 
#'      \item  stk.age.min: Minimum age of the stock (number).
#'      \item  stk.age.max: Maximum age of the stock (number).}
#' Optionals:
#'\itemize{
#'      \item  stk_wt.flq: Mean weight at age of an individual (FLQuant). Required if fl.met.stk_landings.wt.flq is not defined.
#'      \item  stk_n.flq: Numbers at age in the population(FLQuant). Required if Cobb Douglas parameters are not defined.
#'      \item  stk_gB.flq: Biomass growth for the stock modeled in biomass (FLQuant). Required if Cobb Douglas parameters are not defined.}
#'      
#' @return An FLFleetsExt object.
#' 
#
#   Required functions: Create.list.stks.flqa	function, calculate.CDparam
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   Section 1:      FLCatchExt 
#   Section 2:      FLMetierExt 
#   Section 3:      FLFleetExt:
#     3.1               Historical data per fleet 
#     3.2               Projection per fleet 
#     3.3               Historical data per fleet/metier
#     3.4               Projection per fleet/metier
#     3.5               Historical data per fleet/metier/stock
#       3.5.1             Cobb Douglas Parameters: alpha, beta, catch.q
#     3.6               Projection per fleet/metier/stock
#   Section 4:      FLFleetsExt: create fleets
#   Section 5:      Return fleets
#-------------------------------------------------------------------------------


create.fleets.data <- function(yrs,ns,ni,fls.data,stks.data){
  
  fls       <- names(fls.data)
  n.fl      <- length(fls)  
  n.stk     <- length(stks)
  ac        <- as.character
  first.yr <- yrs[["first.yr"]]
  proj.yr  <- yrs[["proj.yr"]]
  last.yr  <- yrs[["last.yr"]]
  
  hist.yrs  <- ac(first.yr:(proj.yr-1))
  proj.yrs  <- ac(proj.yr:last.yr)
  nmy       <- ac(first.yr:last.yr)
  list.stks.unit <- lapply(stks.data, function(ch) grep(pattern="unit", ch, value = TRUE))
  list.stks.age <- lapply(stks.data, function(ch) grep(pattern="age", ch, value = TRUE))
  list.stks.flqa <-  create.list.stks.flqa(stks,yrs,ni,ns,list.stks.unit,list.stks.age)  
  list.stks.flq  <-  create.list.stks.flq(stks,yrs,ni,ns,list.stks.unit)  


  #==============================================================================
  #   Section 1:      FLCatchExt
  #==============================================================================
  
  for( i in 1:n.fl){  #loop fleet
    
    nmfl <- fls[i]
    nmfl.mets  <- get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.mets',sep=''), value = TRUE))
    n.fl.met   <- length(nmfl.mets)
    
    for(j in 1: n.fl.met){   #loop metier
      
      nmfl.met      <- nmfl.mets[j]
      nmfl.met.stks <- get(grep(fls.data[[nmfl]],pattern=paste(nmfl.met,'.stks',sep=''), value = TRUE))
      n.fl.met.stks <- length(nmfl.met.stks)
      list.FLCatchExt <- list()
       for( k in 1:n.fl.met.stks){  #loop stock
        
         nmstk       <- nmfl.met.stks[k]
         stk.flqa    <- list.stks.flqa[[nmstk]]
         stk.flq     <- list.stks.flq[[nmstk]]
         stk.catch   <- FLCatchExt(name = nmstk,landings.n = stk.flqa,discards.n=stk.flqa, 
                                   landings.sel=stk.flqa,discards.sel=stk.flqa, landings.wt=stk.flqa,
                                   discards.wt = stk.flqa,
                                   alpha = stk.flqa, beta = stk.flqa, catch.q = stk.flqa)
         list.FLCatchExt[[k]] <- stk.catch
         #assign(paste(nmstk,'.',nmfl.met,'.catch',sep=''),stk.catch)
      }
    names(list.FLCatchExt) <- nmfl.met.stks
    met.stks.catches <- FLCatchesExt(list.FLCatchExt)
    #met.stks.catches   <- FLCatchesExt(sapply(X=paste(nmfl.met.stks,'.',nmfl.met,'.catch',sep=''),
    #                                          FUN=get, envir=sys.frame(which=-1)))
    #names(met.stks.catches) <- nmfl.met.stks
    assign(paste(nmfl,'.',nmfl.met,'.catch',sep=''),met.stks.catches)
    } 
  }
  
  
  #==============================================================================
  #   Section 2:      FLMetierExt 
  #==============================================================================
  list.FLFleet <- list()
  
  for( i in 1: n.fl){  #loop fleet
 
    nmfl       <- fls[i]
    nmfl.mets  <- get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.mets',sep=''), value = TRUE))
    n.fl.met   <- length(nmfl.mets)   
    efs        <- 0
    list.FLMetierExt <- list()
      for( j in 1:n.fl.met){ #loop metier
        
        nmfl.met<- nmfl.mets[j]
        
        fl.met  <- FLMetierExt(name = nmfl.met, 
                          catches = get(paste(nmfl,'.',nmfl.met,'.catch',sep='')), 
                          effshare = FLQuant(dim=c(1,length(nmy),1,ns),iter=ni, dimnames=list(age='all',year=nmy)),
                          vcost = FLQuant(dim=c(1,length(nmy),1,ns),iter=ni, dimnames=list(age='all',year=nmy)))
        
        list.FLMetierExt [[j]]<- fl.met  
        #assign(paste('fl.',nmfl.met,sep=''), fl.met)
      }
 
    names(list.FLMetierExt) <- nmfl.mets
    
   #==============================================================================
   #   Section 3:      Create FLFleetExt 
   #==============================================================================
   
    cat('\n')
    cat('=============', nmfl,'fleet','=============\n')
    
    met.s <- FLMetiersExt(list.FLMetierExt)
     
  #  met.s <- FLMetiersExt(sapply(paste('fl.',nmfl.mets,sep=''), 
  #                                                       FUN=get, envir=sys.frame(which=-1)))    
    fleet <- FLFleetExt(name = nmfl, 
                           metiers = met.s,
                           effort   = FLQuant(dim=c(1,length(nmy),1,ns),iter= ni, dimnames=list(age='all',year=nmy)),
                           fcost    = FLQuant(dim=c(1,length(nmy),1,ns),iter= ni, dimnames=list(age='all',year=nmy)),
                           capacity = FLQuant(dim=c(1,length(nmy),1,ns),iter= ni, dimnames=list(age='all',year=nmy)),
                           crewshare= FLQuant(dim=c(1,length(nmy),1,ns),iter= ni, dimnames=list(age='all',year=nmy)))
                       
     names(fleet@metiers) <- nmfl.mets

     #-----------------------------------------------------------------------------
     #   Section 3.1:      Historical data per fleet 
     #-----------------------------------------------------------------------------
        
      fl.effort    <- get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'_effort.flq',sep=''), value = TRUE)) 
      fl.fcost     <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'_fcost.flq',sep=''), value = TRUE),envir=as.environment(1)) 
      if(length(fl.fcost)==0) fl.fcost  <- NA
      fl.crewshare     <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'_crewshare.flq',sep=''), value = TRUE),envir=as.environment(1)) 
      if(length(fl.crewshare)==0) fl.crewshare  <- NA
      fl.capacity  <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'_capacity.flq',sep=''), value = TRUE),envir=as.environment(1)) 
      if(length(fl.capacity)==0) fl.capacity  <- NA
      
      # Check dimension names and set units to the fleet object
      if(is.FLQuant(fl.effort)){
        log.dim <- equal.flq.Dimnames(lflq=list(fl.effort,fleet@effort[,hist.yrs]),2)
        if(!log.dim)stop('in effort dimension names \n')
        if(!(any(dim(fl.effort)[4]==c(1,ns))))stop('in effort number of seasons 1 or ns')
        if(!(any(dim(fl.effort)[6]==c(1,ni))))stop('in effort number of iterations 1 or ni')
        units(fleet)$effort <- units(fl.effort)}
      
        
      if(!is.na(fl.fcost)){
        fl.fcost <- fl.fcost[[1]]
        log.dim <- equal.flq.Dimnames(lflq=list(fl.fcost,fleet@fcost[,hist.yrs]),2)
        if(!log.dim)stop('in fcost dimension names  \n')
        if(!(any(dim(fl.fcost)[4]==c(1,ns))))stop('in fcost number of seasons 1 or ns')
        if(!(any(dim(fl.fcost)[6]==c(1,ni))))stop('in fcost number of iterations 1 or ni')
        units(fleet)$fcost <- units(fl.fcost)}
      
        
      if(!is.na(fl.capacity)){
        fl.capacity <- fl.capacity[[1]]
        log.dim <- equal.flq.Dimnames(lflq=list(fl.capacity,fleet@capacity[,hist.yrs]),2)
        if(!log.dim)stop('in capacity dimension names  \n')
        if(!(any(dim(fl.capacity)[4]==c(1,ns))))stop('in capacity number of seasons 1 or ns')
        if(!(any(dim(fl.capacity)[6]==c(1,ni))))stop('in capacity number of iterations 1 or ni')
        units(fleet)$capacity <- units(fl.capacity)}
    
        
      if(!is.na(fl.crewshare)){
        fl.crewshare <- fl.crewshare[[1]]
        log.dim <- equal.flq.Dimnames(lflq=list(fl.crewshare,fleet@crewshare[,hist.yrs]),2)
        if(!log.dim)stop('in crewshare dimension names \n')
        if(!(any(dim(fl.crewshare)[4]==c(1,ns))))stop('in crewshare number of seasons 1 or ns')
        if(!(any(dim(fl.crewshare)[6]==c(1,ni))))stop('in crewshare number of iterations 1 or ni')
        units(fleet)$crewshare <- units(fl.crewshare)}

  
    effort(fleet)[,hist.yrs]      <-  fl.effort 
    fleet@fcost[,hist.yrs]        <-  fl.fcost
    fleet@capacity[,hist.yrs]     <-  fl.capacity
    fleet@crewshare[,hist.yrs]    <-  fl.crewshare
    
     #-----------------------------------------------------------------------------
     #   Section 3.2:      Projection per fleet 
     #-----------------------------------------------------------------------------
    fl.proj.avg.yrs    <- ac(get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'_proj.avg.yrs',sep=''), value = TRUE)))

    effort(fleet)[,proj.yrs,]   <-  yearMeans(effort(fleet)[,fl.proj.avg.yrs,])
    fleet@fcost[,proj.yrs,]     <-  yearMeans(fleet@fcost[,fl.proj.avg.yrs,])
    fleet@capacity[,proj.yrs,]  <-  yearMeans(fleet@capacity[,fl.proj.avg.yrs,])
    fleet@crewshare[,proj.yrs,] <-  yearMeans(fleet@crewshare[,fl.proj.avg.yrs,])
    
    if(any(is.na(effort(fleet)[,fl.proj.avg.yrs]))) { 
      cat('warning: NA-s in effort for average years \n')
      if(any(is.na(effort(fleet)[,proj.yrs])))
        cat('warning: all NA-s in effort for projection years \n')
    }
    
    if(any(is.na(fleet@fcost[,fl.proj.avg.yrs]))) { 
      cat('warning: NA-s in fcost for average years \n')
      if(any(is.na(fleet@fcost[,proj.yrs])))
        cat('warning: all NA-s in fcost for projection years \n')
    }
    
    if(any(is.na(fleet@capacity[,fl.proj.avg.yrs]))) {
      cat('warning: NA-s in capacity for average years \n')
      if(any(is.na(fleet@capacity[,proj.yrs])))
        cat('warning: all NA-s in capacity for projection years \n')
    }
    
    if(any(is.na(fleet@crewshare[,fl.proj.avg.yrs]))) { 
      cat('warning: NA-s in crewshare for average years \n')
      if (any(is.na(fleet@crewshare[,proj.yrs])))
        cat('warning: all NA-s in crewshare for projection years \n')
    }

    all.efs <- numeric(ns) # counter for adding all the effortshares by metier of one fleet for each season
    all.efs[] <- 0
      for (j in 1:n.fl.met){ #loop metier
      
       
        nmfl.met       <- nmfl.mets[j]
        cat('---------------------', nmfl,'fleet,',nmfl.met,'------------------\n')
        
        #-----------------------------------------------------------------------------
        #   3.3     Historic data per fleet/metier
        #-----------------------------------------------------------------------------
       fl.met.effshare    <- get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'_effshare.flq',sep=''), value = TRUE)) 
       # if(length(fl.met.effshare)==0) fl.met.effshare  <- NA
       fl.met.vcost       <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'_vcost.flq',sep=''), value = TRUE),envir=as.environment(1)) 
       if(length(fl.met.vcost)==0) fl.met.vcost  <- NA
       
        # Check dimension names in years and set the units
        if(is.FLQuant(fl.met.effshare)){
          log.dim <- equal.flq.Dimnames(lflq=list(fl.met.effshare, fleet@metiers[[nmfl.met]]@effshare[,hist.yrs]),2)
          if(!log.dim)stop('in effshare dimensions  \n')
          if(!(any(dim(fl.met.effshare)[4]==c(1,ns))))stop('in effshare number of seasons 1 or ns')
          if(!(any(dim(fl.met.effshare)[6]==c(1,ni))))stop('in effshare number of iterations 1 or ni')
          units(fleet@metiers[[nmfl.met]])$effshare <- units(fl.met.effshare)}
        
        if(!is.na(fl.met.vcost)){
          fl.met.vcost <- fl.met.vcost[[1]]
          log.dim <- equal.flq.Dimnames(lflq=list(fl.met.vcost,fleet@metiers[[nmfl.met]]@vcost[,hist.yrs]),2)
          if(!log.dim)stop('in vcost dimensions  \n')
          if(!(any(dim(fl.met.vcost)[4]==c(1,ns))))stop('in vcost number of seasons 1 or ns')
          if(!(any(dim(fl.met.vcost)[6]==c(1,ni))))stop('in vcost number of iterations 1 or ni')
          units(fleet@metiers[[nmfl.met]])$vcost <- units(fl.met.vcost)}
        
        fleet@metiers[[nmfl.met]]@effshare[,hist.yrs] <- fl.met.effshare
        fleet@metiers[[nmfl.met]]@vcost[,hist.yrs]    <- fl.met.vcost
        
        #-----------------------------------------------------------------------------
        #   3.4     Projection per fleet/metier
        #-----------------------------------------------------------------------------

        fl.met.proj.avg.yrs <- ac(get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'_proj.avg.yrs',sep=''), value = TRUE)))

        #   projection effshare
        # NOTE: the sum of all effshare must be one

        for(ss in 1:ns){
          met.sefs <- yearMeans(fl.met.effshare[,fl.met.proj.avg.yrs,,ss])
          if (j==n.fl.met & n.fl.met>1){
            fleet@metiers[[nmfl.met]]@effshare[, proj.yrs,, ss] <- 1 - all.efs[ss]
          } else {
            fleet@metiers[[nmfl.met]]@effshare[,proj.yrs,,ss] <-  met.sefs
          } 
          all.efs[ss] <- all.efs[ss] + met.sefs
        }
        
        if(any(is.na(fl.met.effshare[,fl.met.proj.avg.yrs]))){
          cat('warning: NA-s in effshare for average years \n')
          if(any(is.na(fl.met.effshare[,proj.yrs])))
            cat('warning: all NA-s in effshare for projection years \n')
        }
        
        #   projection vcost
        
        fleet@metiers[[nmfl.met]]@vcost[,proj.yrs,]  <- yearMeans(fleet@metiers[[nmfl.met]]@vcost[,fl.met.proj.avg.yrs,])
        if(any(is.na(fleet@metiers[[nmfl.met]]@vcost[,fl.met.proj.avg.yrs]))) {
          cat('warning: NA-s in vcost for average years \n')
          if(any(is.na(fleet@metiers[[nmfl.met]]@vcost[,proj.yrs])))
            cat('warning: all NA-s in vcost for projection years \n')
        }
        
        
        nmfl.met.stks <- get(paste(nmfl,'.',nmfl.met,'.stks',sep=''))
        n.fl.met.stks <- length(nmfl.met.stks)
        
        for(k in 1:n.fl.met.stks){
    
          nmfl.met.stk <- nmfl.met.stks[k]
          stk.unit     <- get(paste(nmfl.met.stk,'.unit',sep=""))
          flqa.stk     <- list.stks.flqa[[nmfl.met.stk]]
          flq.stk      <- list.stks.flq[[nmfl.met.stk]]
          
          cat('---------------------', nmfl,'fleet,',nmfl.met,' metier,',nmfl.met.stk,'stock','---------------------\n')

          #-----------------------------------------------------------------------------
          #   3.5     Historic data per fleet/metier/stock
          #-----------------------------------------------------------------------------
          landings.n    <- get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_landings.n.flq',sep=''), value = TRUE)) 
          landings.wt   <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_landings.wt.flq',sep=''), value = TRUE),envir=as.environment(1)) 
          if(length(landings.wt)==0) landings.wt <- NA
          discards.wt   <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_discards.wt.flq',sep=''), value = TRUE),envir=as.environment(1)) 
          if(length(discards.wt)==0) discards.wt <- NA
          discards.n   <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_discards.n.flq',sep=''), value = TRUE),envir=as.environment(1)) 
          if(length(discards.n)==0) discards.n <- NA
          price   <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_price.flq',sep=''), value = TRUE),envir=as.environment(1)) 
          if(length(price)==0) price <- NA

          #create the variables with the right dimensions
          
          fl.met.stk.landings.n   <- flqa.stk
          fl.met.stk.discards.n   <- flqa.stk
          fl.met.stk.landings.wt  <- flqa.stk
          fl.met.stk.discards.wt  <- flqa.stk
          fl.met.stk.landings     <- flq.stk
          fl.met.stk.discards     <- flq.stk
          fl.met.stk.landings.sel <- flqa.stk
          fl.met.stk.discards.sel <- flqa.stk
          fl.met.stk.price   <- flqa.stk
          fl.met.stk.alpha   <- flqa.stk
          fl.met.stk.beta    <- flqa.stk
          fl.met.stk.catch.q <- flqa.stk
          
         
          if(all(is.na(landings.n))){
            stop('all NA-s in landings.n')} 
          log.dim <- equal.flq.Dimnames(lflq=list(landings.n,fl.met.stk.landings.n[,hist.yrs]),1:2)
          if(!log.dim)stop('in landings.n dimensions')
          if(!(any(dim(landings.n)[3]==c(1,stk.unit))))stop('in landings.n number of stock units 1 or stk.unit')
          if(!(any(dim(landings.n)[4]==c(1,ns))))stop('in landings.n number of seasons 1 or ns')
          if(!(any(dim(landings.n)[6]==c(1,ni))))stop('in landings.n number of iterations 1 or ni')
          # units(fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]])$landings.n <- units(landings.n)
          units(fl.met.stk.landings.n) <- units(landings.n)
          
          if(!all(is.na(landings.wt))){
            landings.wt <- landings.wt[[1]] 
          log.dim <- equal.flq.Dimnames(lflq=list(landings.wt,fl.met.stk.landings.wt[,hist.yrs]),1:2)
          if(!log.dim)stop('in landings.wt dimensions')
          if(!(any(dim(landings.wt)[3]==c(1,stk.unit))))stop('in landings.wt number of stock units 1 or stk.unit')
          if(!(any(dim(landings.wt)[4]==c(1,ns))))stop('in landings.wt number of seasons 1 or ns')
          if(!(any(dim(landings.wt)[6]==c(1,ni))))stop('in landings.wt number of iterations 1 or ni')
          }else{    
            landings.wt <- get(grep(stks.data[[nmfl.met.stk]],pattern=paste(nmfl.met.stk,'_wt.flq',sep=''), value = TRUE)) 
          } 
          # units(fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]])$landings.wt <- units(landings.wt)
          units(fl.met.stk.landings.wt) <- units(landings.wt)
          
          if(!all(is.na(discards.n))){
            discards.n <- discards.n[[1]]
            log.dim <- equal.flq.Dimnames(lflq=list(discards.n,fl.met.stk.discards.n[,hist.yrs]),1:2)
            if(!log.dim)stop('in discards.n dimensions \n')
            fl.met.stk.discards.n [,hist.yrs] <- discards.n
            if(!(any(dim(discards.n)[3]==c(1,stk.unit))))stop('in discards.n number of stock units 1 or stk.unit')
            if(!(any(dim(discards.n)[4]==c(1,ns))))stop('in discards.n number of seasons 1 or ns')
            if(!(any(dim(discards.n)[6]==c(1,ni))))stop('in discards.n number of iterations 1 or ni')
          }else{    
            cat('warning: all NA-s in discards.n \n')   
          } 
          # units(fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]])$discards.n <- units(discards.n)
          units(fl.met.stk.discards.n) <- 't'# units(discards.n) # dor: there was a bug and I've put 't' just to make it running but it should be corrected
          
          if(!all(is.na(discards.wt))){
            discards.wt <- discards.wt[[1]]
            log.dim <- equal.flq.Dimnames(lflq=list(discards.wt,fl.met.stk.discards.wt[,hist.yrs]),1:2)
            if(!log.dim)stop('in discards.wt dimensions \n')
            fl.met.stk.discards.wt [,hist.yrs] <- discards.wt
            if(!(any(dim(discards.wt)[3]==c(1,stk.unit))))stop('in discards.wt number of stock units 1 or stk.unit')
            if(!(any(dim(discards.wt)[4]==c(1,ns))))stop('in discards.wt number of seasons 1 or ns')
            if(!(any(dim(discards.wt)[6]==c(1,ni))))stop('in discards.wt number of iterations 1 or ni')
          }else{    
            discards.wt <- get(grep(stks.data[[nmfl.met.stk]],pattern=paste(nmfl.met.stk,'_wt.flq',sep=''), value = TRUE))  
          } 
          # units(fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]])$discards.wt <- units(discards.wt)
          units(fl.met.stk.discards.wt) <- 'kg'#units(discards.wt)
          
          
          if(!all(is.na(price))){
            price <- price[[1]]
            log.dim <- equal.flq.Dimnames(lflq=list(price,fl.met.stk.price[,hist.yrs]),1:2)
            if(!log.dim)stop('in price dimensions \n')
            fl.met.stk.price [,hist.yrs] <- price
            if(!(any(dim(price)[3]==c(1,stk.unit))))stop('in price number of stock units 1 or stk.unit')
            if(!(any(dim(price)[4]==c(1,ns))))stop('in price number of seasons 1 or ns')
            if(!(any(dim(price)[6]==c(1,ni))))stop('in price number of iterations 1 or ni')
            # units(fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]])$price <- units(price)
            units(fl.met.stk.price) <- units(price)}

                  
          fl.met.stk.landings.n[,hist.yrs]   <- landings.n
          fl.met.stk.discards.n[,hist.yrs]   <- discards.n
          fl.met.stk.landings.wt[,hist.yrs]  <- landings.wt
          fl.met.stk.discards.wt[,hist.yrs]  <- discards.wt
          fl.met.stk.price[,hist.yrs]  <- price          
          
          #Transformation of NA in landings.n and discards.n in 0
          
          if(all(is.na(fl.met.stk.landings.wt[,hist.yrs]))){
            cat('warning: all NA-s in landings.wt for historic years and will be replaced by 0. \n')
            if(!(any(dim(fl.met.stk.landings.wt)[3]==c(1,stk.unit))))stop('in stk.wt number of stock units 1 or stk.unit')
            if(!(any(dim(fl.met.stk.landings.wt)[4]==c(1,ns))))stop('in stk.wt number of seasons 1 or ns')
            if(!(any(dim(fl.met.stk.landings.wt)[6]==c(1,ni))))stop('in stk.wt number of iterations 1 or ni')}
          
          fl.met.stk.discards.n[,hist.yrs][is.na(fl.met.stk.discards.n[,hist.yrs])] <- 0 
          fl.met.stk.discards.wt[,hist.yrs][is.na(fl.met.stk.discards.wt[,hist.yrs])] <- 0 
          fl.met.stk.landings.n[,hist.yrs][is.na(fl.met.stk.landings.n[,hist.yrs])] <- 0  
          fl.met.stk.landings.wt[,hist.yrs][is.na(fl.met.stk.landings.wt[,hist.yrs])] <- 0  
          
          fl.met.stk.discards     <- unitSums(quantSums(fl.met.stk.discards.n*fl.met.stk.discards.wt))
          fl.met.stk.landings     <- unitSums(quantSums(fl.met.stk.landings.n*fl.met.stk.landings.wt))
          fl.met.stk.landings.sel <- fl.met.stk.landings.n/(fl.met.stk.landings.n+fl.met.stk.discards.n)
          fl.met.stk.discards.sel <- fl.met.stk.discards.n/(fl.met.stk.landings.n+fl.met.stk.discards.n)
            
          if(any((fl.met.stk.landings.n+fl.met.stk.discards.n)==0, na.rm=TRUE)){
            fl.met.stk.landings.sel[(fl.met.stk.landings.n+fl.met.stk.discards.n)==0] <- 1
            fl.met.stk.discards.sel[(fl.met.stk.landings.n+fl.met.stk.discards.n)==0] <- 0}
          
          #-----------------------------------------------------------------------------
          #   3.5.1     Cobb Douglas Parameters: alpha, beta, catch.q
          #-----------------------------------------------------------------------------

          alpha       <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_alpha.flq',sep=''), value = TRUE),envir=as.environment(1))                  
          if(length(alpha)==0) {alpha <- NA
          }else{alpha <- alpha[[1]]}
          beta      <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_beta.flq',sep=''), value = TRUE),envir=as.environment(1))                  
          if(length(beta)==0) {beta <- NA
          }else{beta <- beta[[1]]}
          catch.q     <- mget(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_catch.q.flq',sep=''), value = TRUE),envir=as.environment(1))                  
          if(length(catch.q)==0) {catch.q <- NA
          }else{catch.q <- catch.q[[1]]}
       
          fl.met.stk.alpha[]     <- alpha
          fl.met.stk.beta[]        <- beta
          fl.met.stk.catch.q[]     <- catch.q
          
          if(all(is.na(alpha)) || all(is.na(beta)) || all(is.na(catch.q))){
            CDpar.calc <- TRUE
            stk.n       <- get(grep(stks.data[[ nmfl.met.stk]],pattern=paste(nmfl.met.stk,'_n.flq',sep=''), value = TRUE)) 
            stk.n[,hist.yrs][is.na(stk.n[,hist.yrs])] <- 0
            fl.effort[,hist.yrs][is.na(fl.effort[,hist.yrs])] <- 0
            stk.age.min <- get(grep(stks.data[[ nmfl.met.stk]],pattern=paste(nmfl.met.stk,'.age.min',sep=''), value = TRUE)) 
            stk.age.max <- get(grep(stks.data[[ nmfl.met.stk]],pattern=paste(nmfl.met.stk,'.age.max',sep=''), value = TRUE)) 
            
            largs <- NULL
            
            if(all(is.na(c(stk.age.max,stk.age.min)))){
              largs$stk.gB <- get(grep(stks.data[[ nmfl.met.stk]],pattern=paste(nmfl.met.stk,'_gB.flq',sep=''), value = TRUE)) 
            }else{
              if(length(stk.age.min:stk.age.max)==1 ){
              largs$stk.gB <- get(grep(stks.data[[ nmfl.met.stk]],pattern=paste(nmfl.met.stk,'_gB.flq',sep=''), value = TRUE)) 
              }else{ 
              largs$stk.m <- get(grep(stks.data[[ nmfl.met.stk]],pattern=paste(nmfl.met.stk,'_m.flq',sep=''), value = TRUE))
              }}
            
            CD_param <- calculate.CDparam(stk.n, fl.met.stk.landings.n,fl.met.stk.discards.n,
                                         fl.effort,fl.met.effshare,stk.age.min,stk.age.max,
                                         flqa.stk,flq.stk,largs)
          
            fl.met.stk.alpha    <- CD_param[['alpha']]
            fl.met.stk.beta     <- CD_param[['beta']]
            fl.met.stk.catch.q  <- CD_param[['catch.q']]
            
            fl.met.stk.proj.avg.yrs <- ac(get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_proj.avg.yrs',sep=''), value = TRUE)))                  
            
            fl.met.stk.alpha[,proj.yrs,]   <- yearMeans(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs,])
            fl.met.stk.beta[,proj.yrs,]    <- yearMeans(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs,])
            fl.met.stk.catch.q[,proj.yrs,] <- yearMeans(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs,])
            
            }else{
              CDpar.calc <- FALSE
              # Check dimension names
              log.dim <- equal.flq.Dimnames(lflq=list(alpha,beta,catch.q,flqa.stk),1:2)
              if(!log.dim)stop('in alpha,beta or catch.q dimensions \n')
              if(!(any(dim(fl.met.stk.alpha)[3]==c(1,stk.unit))))stop('in alpha number of stock units 1 or stk.unit')
              if(!(any(dim(fl.met.stk.alpha)[4]==c(1,ns))))stop('in alpha number of seasons 1 or ns')
              if(!(any(dim(fl.met.stk.alpha)[6]==c(1,ni))))stop('in alpha number of iterations 1 or ni')
              if(!(any(dim(fl.met.stk.beta)[3]==c(1,stk.unit))))stop('in beta number of stock units 1 or stk.unit')
              if(!(any(dim(fl.met.stk.beta)[4]==c(1,ns))))stop('in beta number of seasons 1 or ns')
              if(!(any(dim(fl.met.stk.beta)[6]==c(1,ni))))stop('in beta number of iterations 1 or ni')       
              if(!(any(dim(fl.met.stk.catch.q)[3]==c(1,stk.unit))))stop('in catch.q number of stock units 1 or stk.unit')
              if(!(any(dim(fl.met.stk.catch.q)[4]==c(1,ns))))stop('in catch.q number of seasons 1 or ns')
              if(!(any(dim(fl.met.stk.catch.q)[6]==c(1,ni))))stop('in catch.q number of iterations 1 or ni') 
              
              # if catch.q historical values are set --> values are required also for projection years
              if (any(is.na(fl.met.stk.catch.q[, proj.yrs]))) {
                cat('NA-s in catch.q projection. As historical values were set, 
                    then projection values must also be set as they are not estimated. \n')
              }
            }
          
          #-----------------------------------------------------------------------------
          #   3.6     Projection per fleet/metier/stock
          #-----------------------------------------------------------------------------
          
          fl.met.stk.proj.avg.yrs <- ac(get(grep(fls.data[[nmfl]],pattern=paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_proj.avg.yrs',sep=''), value = TRUE)))                  

          if(any(is.na(fl.met.stk.landings.sel[, fl.met.stk.proj.avg.yrs]))) {
            cat('warning: NA-s in landings.sel for average years and will be replaced by 1. \n')
            fl.met.stk.landings.sel[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.landings.sel[, fl.met.stk.proj.avg.yrs])] <- 1
            # if (any(is.na(fl.met.stk.landings.sel[, proj.yrs])))
            #   cat('warning: all NA-s in landings.sel projection. \n')
          }
          fl.met.stk.landings.sel[,proj.yrs,] <- yearMeans(fl.met.stk.landings.sel[,fl.met.stk.proj.avg.yrs,])
          
          
          if(any(is.na(fl.met.stk.discards.sel[, fl.met.stk.proj.avg.yrs]))) {
            cat('warning: NA-s in discards.sel for average years and will be replaced by 0. n')
            fl.met.stk.discards.sel[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.discards.sel[, fl.met.stk.proj.avg.yrs])] <- 0
            # if (any(is.na(fl.met.stk.discards.sel[, proj.yrs])))
            #   cat('warning: all NA-s in discards.sel projection\n')
          }
          fl.met.stk.discards.sel[,proj.yrs,] <- yearMeans(fl.met.stk.discards.sel[,fl.met.stk.proj.avg.yrs,])
          
          fl.met.stk.landings.wt[,proj.yrs,] <- yearMeans(fl.met.stk.landings.wt[, fl.met.stk.proj.avg.yrs,])
          
          fl.met.stk.discards.wt[,proj.yrs,] <- yearMeans(fl.met.stk.discards.wt[, fl.met.stk.proj.avg.yrs,])
          
          fl.met.stk.price[,proj.yrs,]       <- yearMeans(fl.met.stk.price[, fl.met.stk.proj.avg.yrs,])
          
          if (any(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs]<0, na.rm = TRUE)) { 
            stop('Negative values in alpha projection. \n')
          } else if (any(is.na(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs]))) {
            cat('warning: NA-s in alpha for average years and will be replaced by 1. \n')
            fl.met.stk.alpha[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs])] <- 1
            if (CDpar.calc == TRUE) fl.met.stk.alpha[,proj.yrs,] <- yearMeans(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs,])
            # if (any(is.na(fl.met.stk.alpha[, proj.yrs])))
            #   cat('warning: all NA-s in alpha projection \n')
          }
          
          if (any(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs]<0, na.rm = TRUE)) { 
            stop('Negative values in beta projection. \n')
          } else if (any(is.na(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs]))) {
            cat('warning: NA-s in beta for average years and will be replaced by 1. \n')
            fl.met.stk.beta[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs])] <- 1
            if (CDpar.calc == TRUE) fl.met.stk.beta[,proj.yrs,] <- yearMeans(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs,])
            # if (any(is.na(fl.met.stk.beta[, proj.yrs])))
            #   cat('warning: all NA-s in beta projection \n')
          }
          
          if (any(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs]<0, na.rm = TRUE)) { 
            stop('Negative values in catch.q projection. \n')
          } else if (any(is.na(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs]))) {
            cat('warning: NA-s in catch.q for average years \n') 
            if (any(is.na(fl.met.stk.catch.q[, proj.yrs]))) {
              cat('warning: all NA-s in catch.q projection. \n')
            }
          }
          
          landings.n(fleet, metier = nmfl.met, catch = nmfl.met.stk)   <-  fl.met.stk.landings.n
          discards.n(fleet, metier = nmfl.met, catch = nmfl.met.stk)   <-  fl.met.stk.discards.n
          landings(fleet, metier = nmfl.met, catch = nmfl.met.stk)     <-  fl.met.stk.landings
          discards(fleet, metier = nmfl.met, catch = nmfl.met.stk)     <-  fl.met.stk.discards
          landings.sel(fleet, metier = nmfl.met, catch = nmfl.met.stk) <-  fl.met.stk.landings.sel
          discards.sel(fleet, metier = nmfl.met, catch = nmfl.met.stk) <-  fl.met.stk.discards.sel
          landings.wt(fleet, metier = nmfl.met, catch = nmfl.met.stk)  <-  fl.met.stk.landings.wt
          discards.wt(fleet, metier = nmfl.met, catch = nmfl.met.stk)  <-  fl.met.stk.discards.wt
          price(fleet, metier = nmfl.met, catch = nmfl.met.stk)        <-  fl.met.stk.price
          fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]]@alpha      <-  fl.met.stk.alpha
          fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]]@beta       <-  fl.met.stk.beta
          fleet@metiers[[nmfl.met]]@catches[[nmfl.met.stk]]@catch.q    <-  fl.met.stk.catch.q

        }  # loop stock
        
      } # loop metier
    
    # assign(paste(nmfl,'.fleet',sep=''), fleet)
    list.FLFleet[[i]]<- fleet
    
    # #Checking that the sum is close to 1.
    # sum.efsh <- 0
    # sum.yr <- length(proj.yrs)
    # for (ss in 1:ns) {
    #   for (j in 1:n.fl.met) {
    #     nmfl.met <- nmfl.mets[j]
    #     sum.efsh<- sum.efsh+ fleet@metiers[[nmfl.met]]@effshare[, proj.yrs[sum.yr], , ss]  #/all.efs
    #   }
    #   if(abs(sum.efsh-1)>=10^(-3)){ 
    #     stop(paste("The total sum of effshare is not one in season ",ss," and fleet ", fls[i],sep=""))
    #   }
    # }
  } # loop fleet
  
  #==============================================================================
  #   Section 4:     FLFleetsExt: create fleets
  #==============================================================================
  names(list.FLFleet) <- fls
  fleets <- FLFleetsExt(list.FLFleet)
  #fleets        <- FLFleetsExt(sapply(paste(fls,'.fleet',sep=''),FUN=get, envir=sys.frame(which=-1)))
  #names(fleets) <- fls
  
  #==============================================================================
  #   Section 5:           Return
  #==============================================================================
  
  return(fleets)
  
}
