###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.fleets.data
# NOTE #1:      Return FLFleets object called fleets
###############################################################################
#-------------------------------------------------------------------------
#  inputs: 
#
#   (required)
#   first.yr: First year of simulation (number)
#   proj.yr:  First year of projection (number)
#   last.yr:  Last year of projection (number)
#   ni:       Number of iterations (number)
#   ns:	      Number of seasons (number)
#   stks:     Name of all the stocks (vector)
#   stk.unit	Number of units of the stock (number) 
#   stk.min.age:  Minimum age of the stock (number)
#   stk.max.age:  Maximum age of the stock (number)
#   stk_wt.flq:   Weight age age (FLQuant)
#   fls:	        Name of all the fleets (vector)
#   fl.met.stks:	Name of the stocks in the metier 'met' and fleet 'fl' (vector)
#   fl_effort.flq: 'fl' fleet's effort (FLQuant)
#   fl.met_effshare.flq:        'fl' fleet and 'met' metier's effort share (FLQuant)
#   fl.met.stk_landings.n.flq:  'fl' fleet,'met' m?tier and 'stk' stocks landings at age
#   fl_proj.avg.yrs  vector:	  historic years to calculate the average of effort,fcost,crewshare,capacity in 'fl' fleet (vector)
#
#   (optionals)
#   fl_capacity.flq:	  'fl' fleet's capacity (FLQuant)
#   fl_fcost.flq:       'fl' fleet's fixed cost (FLQuant)
#   fl_crewshare.flq:   'fl' fleet's crewshare (FLQuant)
#   fl.met_vcost.flq:	  'fl' fleet and 'met' metier's variable cost (FLQuant)
#   fl.met.stk_discards.n.flq: 'fl' fleet,'met' m?tier and 'stk' stocks discards at age (FLQuant)
#   fl.met.stk_price.flq:	'fl' fleet,'met' m?tier and 'stk' stocks price at age (FLQuant)
#   fl.met.stk_alpha.flq:	'fl' fleet,'met' m?tier and 'stk' stocks Cobb Douglass alpha parameter (FLQuant)
#   fl.met.stk_beta.flq:	'fl' fleet,'met' m?tier and 'stk' stocks Cobb Douglass beta parameter (FLQuant)
#   fl.met.stk_catch.q.flq:  'fl' fleet,'met' m?tier and 'stk' stocks Cobb Douglass catch.q parameter (FLQuant)
#   stk_n.flq:          Abundance at age (if alpha, beta or catch.q is not defined required) (FLQuant) 
#   fl.met_proj.avg.yrs:	    historic years to calculate the average of effshare,vcost for 'fl' fleet and 'met' metier(vector)
#   fl.met.stk_proj.avg.yrs:	historic years to calculate the average of landings.wt, discards.wt,landings.sel,discards.sel
#                                   alpha,beta,catch.q for 'fl' fleet, 'met' metier and 'stk' stock(vector)
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


create.fleets.data <- function(){
  
  n.fl      <- length(fls)  
  n.stk     <- length(stks)
  ac        <- as.character
  hist.yrs  <- ac(first.yr:(proj.yr-1))
  proj.yrs  <- ac(proj.yr:last.yr)
  nmy       <- ac(first.yr:last.yr)
  
  list.stks.flqa <- create.list.stks.flqa ()
  list.stks.flq  <- create.list.stks.flq ()


  #==============================================================================
  #   Section 1:      FLCatchExt
  #==============================================================================
  
  for( i in 1:n.fl){  #loop fleet
    
    nmfl <- fls[i]
    nmfl.mets  <- get(paste(nmfl,'.mets',sep=''))
    n.fl.met   <- length(nmfl.mets)
    
    for(j in 1: n.fl.met){   #loop metier
      
      nmfl.met      <- nmfl.mets[j]
      nmfl.met.stks <- get(paste(nmfl,'.',nmfl.met,'.stks',sep=''))
      n.fl.met.stks <- length(nmfl.met.stks)
      
       for( k in 1:n.fl.met.stks){  #loop stock
        
         nmstk       <- nmfl.met.stks[k]
         stk.flqa    <- list.stks.flqa[[nmstk]]
         stk.flq     <- list.stks.flq[[nmstk]]
         stk.catch   <- FLCatchExt(name = nmstk,landings.n = stk.flqa,discards.n=stk.flqa, 
                                   landings.sel=stk.flqa,discards.sel=stk.flqa, landings.wt=stk.flqa,
                                   discards.wt = stk.flqa,
                                   alpha = stk.flqa, beta = stk.flqa, catch.q = stk.flqa)
         
         assign(paste(nmstk,'.',nmfl.met,'.catch',sep=''),stk.catch)
      }

    met.stks.catches   <- FLCatchesExt(sapply(X=paste(nmfl.met.stks,'.',nmfl.met,'.catch',sep=''),
                                              FUN=get, envir=sys.frame(which=-1)))
    names(met.stks.catches) <- nmfl.met.stks
    assign(paste(nmfl,'.',nmfl.met,'.catch',sep=''),met.stks.catches)
    } 
  }
  
  
  #==============================================================================
  #   Section 2:      FLMetierExt 
  #==============================================================================
  
  for( i in 1: n.fl){  #loop fleet
 
    nmfl       <- fls[i]
    nmfl.mets  <- get(paste(nmfl,'.mets',sep=''))
    n.fl.met   <- length(nmfl.mets)   
    efs        <- 0
    
      for( j in 1:n.fl.met){ #loop metier
        
        nmfl.met<- nmfl.mets[j]
        
        fl.met  <- FLMetierExt(name = nmfl.met, 
                          catches = get(paste(nmfl,'.',nmfl.met,'.catch',sep='')), 
                          effshare = FLQuant(dim=c(1,length(nmy),1,ns),iter=ni, dimnames=list(age='all',year=nmy)),
                          vcost = FLQuant(dim=c(1,length(nmy),1,ns),iter=ni, dimnames=list(age='all',year=nmy)))
          
        assign(paste('fl.',nmfl.met,sep=''), fl.met)
      }

   #==============================================================================
   #   Section 3:      Create FLFleetExt 
   #==============================================================================
   
    cat('\n')
    cat('=============', nmfl,'fleet','=============\n')
     
    met.s <- FLMetiersExt(sapply(paste('fl.',nmfl.mets,sep=''), 
                                                         FUN=get, envir=sys.frame(which=-1)))    
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
        
      fl.effort    <- get(paste(nmfl,'_effort.flq',sep=''))
      fl.fcost     <- mget(paste(nmfl,'_fcost.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
      fl.crewshare <- mget(paste(nmfl,'_crewshare.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
      fl.capacity  <- mget(paste(nmfl,'_capacity.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
      
      # Check dimension names
      if(!is.na(fl.effort)){
        log.dim <- equal.flq.Dimnames(lflq=list(fl.effort,fleet@effort[,hist.yrs]),2)
        if(!log.dim)stop('in effort dimension names \n')
        if(!(any(dim(fl.effort)[4]==c(1,ns))))stop('in effort number of seasons 1 or ns')
        if(!(any(dim(fl.effort)[6]==c(1,ni))))stop('in effort number of iterations 1 or ni')}
      
        
      if(!is.na(fl.fcost)){
        log.dim <- equal.flq.Dimnames(lflq=list(fl.fcost,fleet@fcost[,hist.yrs]),2)
        if(!log.dim)stop('in fcost dimension names  \n')
        if(!(any(dim(fl.fcost)[4]==c(1,ns))))stop('in fcost number of seasons 1 or ns')
        if(!(any(dim(fl.fcost)[6]==c(1,ni))))stop('in fcost number of iterations 1 or ni')}
      
        
      if(!is.na(fl.capacity)){
        log.dim <- equal.flq.Dimnames(lflq=list(fl.capacity,fleet@capacity[,hist.yrs]),2)
        if(!log.dim)stop('in capacity dimension names  \n')
        if(!(any(dim(fl.capacity)[4]==c(1,ns))))stop('in capacity number of seasons 1 or ns')
        if(!(any(dim(fl.capacity)[6]==c(1,ni))))stop('in capacity number of iterations 1 or ni')}
    
        
      if(!is.na(fl.crewshare)){
        log.dim <- equal.flq.Dimnames(lflq=list(fl.crewshare,fleet@crewshare[,hist.yrs]),2)
        if(!log.dim)stop('in crewshare dimension names \n')
        if(!(any(dim(fl.crewshare)[4]==c(1,ns))))stop('in crewshare number of seasons 1 or ns')
        if(!(any(dim(fl.crewshare)[6]==c(1,ni))))stop('in crewshare number of iterations 1 or ni')}
    
  
    effort(fleet)[,hist.yrs]      <-  fl.effort 
    fleet@fcost[,hist.yrs]        <-  fl.fcost
    fleet@capacity[,hist.yrs]     <-  fl.capacity
    fleet@crewshare[,hist.yrs]    <-  fl.crewshare
    
     #-----------------------------------------------------------------------------
     #   Section 3.2:      Projection per fleet 
     #-----------------------------------------------------------------------------
        
    fl.proj.avg.yrs   <- mget(paste(nmfl,'_proj.avg.yrs',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
    fl.proj.avg.yrs <- ac(fl.proj.avg.yrs)
    
    for(ss in 1:ns){
      effort(fleet)[,proj.yrs,,ss]   <-  yearMeans(effort(fleet)[,fl.proj.avg.yrs,,ss])
      fleet@fcost[,proj.yrs,,ss]     <-  yearMeans(fleet@fcost[,fl.proj.avg.yrs,,ss])
      fleet@capacity[,proj.yrs,,ss]  <-  yearMeans(fleet@capacity[,fl.proj.avg.yrs,,ss])
      fleet@crewshare[,proj.yrs,,ss] <-  yearMeans(fleet@crewshare[,fl.proj.avg.yrs,,ss])
    }
    
    if(any(is.na(effort(fleet)[,fl.proj.avg.yrs]))) { 
      cat('warning: all NA-s in effort projection \n')}
    if(any(is.na(fleet@fcost[,fl.proj.avg.yrs]))) { 
      cat('warning: all NA-s in fcost projection \n')}
    if(any(is.na(fleet@capacity[,fl.proj.avg.yrs]))) {
      cat('warning: all NA-s in capacity projection \n')}
    if(any(is.na(fleet@crewshare[,fl.proj.avg.yrs]))) { 
      cat('warning: all NA-s in crewshare projection \n')}
    
      for (j in 1:n.fl.met){ #loop metier
      
        nmfl.met       <- nmfl.mets[j]
        cat('---------------------', nmfl,'fleet,',nmfl.met,'------------------\n')
        
        #-----------------------------------------------------------------------------
        #   3.3     Historic data per fleet/metier
        #-----------------------------------------------------------------------------
        
        fl.met.effshare <- mget(paste(nmfl,'.',nmfl.met,'_effshare.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
        fl.met.vcost    <- mget(paste(nmfl,'.',nmfl.met,'_vcost.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]

        # Check dimension names in years
        
        if(!is.na(fl.met.effshare)){
          log.dim <- equal.flq.Dimnames(lflq=list(fl.met.effshare, fleet@metiers[[nmfl.met]]@effshare[,hist.yrs]),2)
          if(!log.dim)stop('in effshare dimensions  \n')
          if(!(any(dim(fl.met.effshare)[4]==c(1,ns))))stop('in effshare number of seasons 1 or ns')
          if(!(any(dim(fl.met.effshare)[6]==c(1,ni))))stop('in effshare number of iterations 1 or ni')}
        
        if(!is.na(fl.met.vcost)){
          log.dim <- equal.flq.Dimnames(lflq=list(fl.met.vcost,fleet@metiers[[nmfl.met]]@vcost[,hist.yrs]),2)
          if(!log.dim)stop('in vcost dimensions  \n')
          if(!(any(dim(fl.met.vcost)[4]==c(1,ns))))stop('in vcost number of seasons 1 or ns')
          if(!(any(dim(fl.met.vcost)[6]==c(1,ni))))stop('in vcost number of iterations 1 or ni')}
        
        fleet@metiers[[nmfl.met]]@effshare[,hist.yrs] <- fl.met.effshare
        fleet@metiers[[nmfl.met]]@vcost[,hist.yrs]    <- fl.met.vcost
        
        #-----------------------------------------------------------------------------
        #   3.4     Projection per fleet/metier
        #-----------------------------------------------------------------------------
        
        fl.met.proj.avg.yrs <- ac(get(paste(nmfl,'.',nmfl.met,'_proj.avg.yrs',sep='')))
        
        #   projection effshare
        # NOTE: the sum of all effshare must be one
        
        all.efs <- 0
        for(ss in 1:ns){
          all.efs <- all.efs+ yearMeans(fl.met.effshare[,fl.met.proj.avg.yrs,,ss])} 

        for(ss in 1:ns){
          if (!all.efs == 0){
          fleet@metiers[[nmfl.met]]@effshare[,proj.yrs,,ss] <-  yearMeans(fl.met.effshare[,fl.met.proj.avg.yrs,,ss])/all.efs 
          }else{
          fleet@metiers[[nmfl.met]]@effshare[,proj.yrs,,ss] <- 0
          }
        }
        
        if(any(is.na(fl.met.effshare[,fl.met.proj.avg.yrs]))){
          stop('warning: NA in effshare projection')}
        
        #   projection vcost
        
        if(all(is.na(fl.met.vcost))){
          cat('warning: NA in vcost projection \n')}
        for(ss in 1:ns){
          fleet@metiers[[nmfl.met]]@vcost[,proj.yrs,,ss]  <- yearMeans(fleet@metiers[[nmfl.met]]@vcost[,fl.met.proj.avg.yrs,,ss])
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
    
          landings.n  <- get(paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_landings.n.flq',sep=''))
          landings.wt <- get(paste(nmfl.met.stk,'_wt.flq',sep='')) 
          discards.wt <- get(paste(nmfl.met.stk,'_wt.flq',sep=''))          
          discards.n  <- mget(paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_discards.n.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
          price       <- mget(paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_price.flq',sep=''),envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
          
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
          
          fl.met.stk.landings.n[,hist.yrs]   <- landings.n
          fl.met.stk.discards.n[,hist.yrs]   <- discards.n
          fl.met.stk.landings.wt[,hist.yrs]  <- landings.wt
          fl.met.stk.discards.wt[,hist.yrs]  <- discards.wt
          fl.met.stk.price[,hist.yrs]  <- price
        
          if(all(is.na(landings.n))){
            stop('all NA-s in landings.n')} 
          
          log.dim <- equal.flq.Dimnames(lflq=list(landings.n,fl.met.stk.landings.n[,hist.yrs]),1:2)
          if(!log.dim)stop('in landings.n dimensions')
          if(!(any(dim(landings.n)[3]==c(1,stk.unit))))stop('in landings.n number of stock units 1 or stk.unit')
          if(!(any(dim(landings.n)[4]==c(1,ns))))stop('in landings.n number of seasons 1 or ns')
          if(!(any(dim(landings.n)[6]==c(1,ni))))stop('in landings.n number of iterations 1 or ni')
        
          
          if(!all(is.na(discards.n))){
            log.dim <- equal.flq.Dimnames(lflq=list(discards.n,fl.met.stk.discards.n[,hist.yrs]),1:2)
            if(!log.dim)stop('in discards.n dimensions \n')
            fl.met.stk.discards.n [,hist.yrs] <- discards.n
            if(!(any(dim(discards.n)[3]==c(1,stk.unit))))stop('in discards.n number of stock units 1 or stk.unit')
            if(!(any(dim(discards.n)[4]==c(1,ns))))stop('in discards.n number of seasons 1 or ns')
            if(!(any(dim(discards.n)[6]==c(1,ni))))stop('in discards.n number of iterations 1 or ni')
          }else{    
            cat('warning: all NA-s in discards.n \n')   
          }  

          if(!all(is.na(price))){
            log.dim <- equal.flq.Dimnames(lflq=list(price,fl.met.stk.price[,hist.yrs]),1:2)
            if(!log.dim)stop('in price dimensions \n')
            fl.met.stk.price [,hist.yrs] <- price
            if(!(any(dim(price)[3]==c(1,stk.unit))))stop('in price number of stock units 1 or stk.unit')
            if(!(any(dim(price)[4]==c(1,ns))))stop('in price number of seasons 1 or ns')
            if(!(any(dim(price)[6]==c(1,ni))))stop('in price number of iterations 1 or ni')}
          
          
          #Transformation of NA in landings.n and discards.n in 0
          
          if(all(is.na(fl.met.stk.landings.wt[,hist.yrs]))){
            stop('warning: all NA-s in landings.wt \n')
            if(!(any(dim(fl.met.stk.landings.wt)[3]==c(1,stk.unit))))stop('in stk.wt number of stock units 1 or stk.unit')
            if(!(any(dim(fl.met.stk.landings.wt)[4]==c(1,ns))))stop('in stk.wt number of seasons 1 or ns')
            if(!(any(dim(fl.met.stk.landings.wt)[6]==c(1,ni))))stop('in stk.wt number of iterations 1 or ni')}
          
          fl.met.stk.discards.n[,hist.yrs][is.na(fl.met.stk.discards.n[,hist.yrs])]   <- 0
          fl.met.stk.landings.n[,hist.yrs][is.na(fl.met.stk.landings.n[,hist.yrs])]   <- 0
          fl.met.stk.landings.wt[,hist.yrs][is.na(fl.met.stk.landings.wt[,hist.yrs])] <- 0
          fl.met.stk.discards.wt[,hist.yrs][is.na(fl.met.stk.discards.wt[,hist.yrs])] <- 0
          
          fl.met.stk.discards     <- unitSums(quantSums(fl.met.stk.discards.n*fl.met.stk.discards.wt))
          fl.met.stk.landings     <- unitSums(quantSums(fl.met.stk.landings.n*fl.met.stk.landings.wt))
          fl.met.stk.landings.sel <- fl.met.stk.landings.n/(fl.met.stk.landings.n+fl.met.stk.discards.n)
          fl.met.stk.discards.sel <- fl.met.stk.discards.n/(fl.met.stk.landings.n+fl.met.stk.discards.n)
            
          if(any((fl.met.stk.landings.n+fl.met.stk.discards.n)==0, na.rm=TRUE)){
            fl.met.stk.landings.sel[(fl.met.stk.landings.n+fl.met.stk.discards.n)==0] <- 0
            fl.met.stk.discards.sel[(fl.met.stk.landings.n+fl.met.stk.discards.n)==0] <- 0}
          
          #-----------------------------------------------------------------------------
          #   3.5.1     Cobb Douglas Parameters: alpha, beta, catch.q
          #-----------------------------------------------------------------------------

          alpha       <- paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_alpha.flq',sep='')
          beta        <- paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_beta.flq',sep='')
          catch.q     <- paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_catch.q.flq',sep='')
        
          fl.met.stk.alpha[]       <- mget(alpha,envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
          fl.met.stk.beta[]        <- mget(beta,envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
          fl.met.stk.catch.q[]     <- mget(catch.q,envir=as.environment(-1),ifnotfound=NA,inherits=TRUE)[[1]]
          
          if(any(!exists(alpha) || !exists(beta) || !exists(catch.q))){
            stk.n       <- get(paste(nmfl.met.stk,'_n.flq',sep=''))
            stk.n[,hist.yrs][is.na(stk.n[,hist.yrs])] <- 0
            fl.effort[,hist.yrs][is.na(fl.effort[,hist.yrs])] <- 0
            stk.age.min <- get(paste(nmfl.met.stk,'.age.min',sep=''))
            stk.age.max <- get(paste(nmfl.met.stk,'.age.max',sep=''))

            CD_param <- calculate.CDparam(stk.n, fl.met.stk.landings.n,fl.met.stk.discards.n,
                                         fl.effort,fl.met.effshare,stk.age.min,stk.age.max,
                                         flqa.stk,flq.stk)
          
            fl.met.stk.alpha    <- CD_param[['alpha']]
            fl.met.stk.beta     <- CD_param[['beta']]
            fl.met.stk.catch.q  <- CD_param[['catch.q']]
            
          }else{
            # Check dimension names
            log.dim <- equal.flq.Dimnames(lflq=list(get(alpha),get(beta),get(catch.q),flqa.stk),1:2)
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
            }
          
          #-----------------------------------------------------------------------------
          #   3.6     Projection per fleet/metier/stock
          #-----------------------------------------------------------------------------
          
          fl.met.stk.proj.avg.yrs <- ac(get(paste(nmfl,'.',nmfl.met,'.',nmfl.met.stk,'_proj.avg.yrs',sep='') ))
          
          if(any(is.na(fl.met.stk.landings.sel[, fl.met.stk.proj.avg.yrs]))){
            cat('warning: all NA-s in landings.sel projection \n')
            fl.met.stk.landings.sel[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.landings.sel[, fl.met.stk.proj.avg.yrs])]} <- 0
            
          if(any(is.na(fl.met.stk.discards.sel[, fl.met.stk.proj.avg.yrs]))){
            cat('warning: all NA-s in discards.sel projection \n')
            fl.met.stk.discards.sel[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.discards.sel[, fl.met.stk.proj.avg.yrs])] <- 0}
            
          if(any(is.na(fl.met.stk.price[, fl.met.stk.proj.avg.yrs]))){
            cat('warning: all NA-s in price projection \n')
            fl.met.stk.price[,fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.price[, fl.met.stk.proj.avg.yrs])] <- 0}

          if(any(is.na(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs]))){
            cat('warning: all NA-s in alpha projection \n')
            fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs])] <- 0
            }else{if(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs]<0)
            stop('<0 values in alpha projection \n')}
             
          if(any(is.na(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs]))){
            cat('warning: all NA-s in beta projection \n')
            fl.met.stk.beta[, fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs])] <- 0
            }else{ if(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs][, fl.met.stk.proj.avg.yrs]<0)
            stop('<0 values in beta projection \n')}
                                
          if(any(is.na(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs]))){
            cat('warning: all NA-s in catch.q projection \n')
            fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs][is.na(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs])] <- 0
            }else{ if(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs]<0)
            stop('<0 values in catch.q projection \n')}
          
          for(ss in 1:ns){
            for(unit in 1:stk.unit){
              fl.met.stk.landings.sel[,proj.yrs,unit,ss] <- yearMeans(fl.met.stk.landings.sel[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.discards.sel[,proj.yrs,unit,ss] <- yearMeans(fl.met.stk.discards.sel[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.landings.wt[,proj.yrs,unit,ss]  <- yearMeans(fl.met.stk.landings.wt[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.discards.wt[,proj.yrs,unit,ss]  <- yearMeans(fl.met.stk.discards.wt[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.price[,proj.yrs,unit,ss]   <- yearMeans(fl.met.stk.price[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.alpha[,proj.yrs,unit,ss]   <- yearMeans(fl.met.stk.alpha[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.beta[,proj.yrs,unit,ss]    <- yearMeans(fl.met.stk.beta[, fl.met.stk.proj.avg.yrs,unit,ss])
              fl.met.stk.catch.q[,proj.yrs,unit,ss] <- yearMeans(fl.met.stk.catch.q[, fl.met.stk.proj.avg.yrs,unit,ss])
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
        
      } #loop metier
    
    assign(paste(nmfl,'.fleet',sep=''), fleet)
    
  }        #loop fleet
  
  #==============================================================================
  #   Section 4:     FLFleetsExt: create fleets
  #==============================================================================
  
  fleets        <- FLFleetsExt(sapply(paste(fls,'.fleet',sep=''),FUN=get, envir=sys.frame(which=-1)))
  names(fleets) <- fls

  #==============================================================================
  #   Section 5:           Return
  #==============================================================================
  
  return(fleets)

}
