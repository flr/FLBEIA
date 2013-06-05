#-------------------------------------------------------------------------------
#                    ** create.observation.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 29/05/2013 11:34:54
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - **  ** : character vector with covariables names
# - **  ** : characted vector with the same length as cvrsnames with the process model followed by each of the covariables. 
#         the first element correspond with the process model of the first covariable in cvrsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCovar' is used for **all** the fleets.    

create.obs.ctrl.t <- function(stksnames, n.stks.inds = NULL, stks.indsnames = NULL, stkObs.models = NULL, indObs.models = NULL, immediate = FALSE,...){

 
    stkObs.models.available <- c("age2ageDat", "age2agePop", "age2bioDat", "age2bioPop",
                                 "bio2bioDat", "bio2bioPop", "NoObsStock", "perfectObs")
  
    indObs.models.available <- c("bioInd", "create.ageInd.ctrl", "NoObsIndex") 
    
    res        <- vector('list', length(stksnames))
    names(res) <- stksnames
    extra.args <- list(...)
    
    if(is.null(n.stks.inds)) n.stks.inds <- rep(0,length(stksnames)) 
    
    if(is.null(stkObs.models)) stkObs.models <- rep('perfectObs', length(stksnames))
    if(is.null(indObs.models)) indObs.models <- rep('NoObsIndex', length(stksnames))
    
    if(length(stkObs.models) < length(stksnames)) stop("'stkObs.models' must be NULL or must have the same length as stksnames'")
    if(length(indObs.models) < length(stksnames)) stop("'indObs.models' must be NULL or must have the same length as stksnames'")
    
    if(!all(stkObs.models %in% stkObs.models.available)){ 
      wmod <- unique(stkObs.models[which(!(stkObs.models %in% stkObs.models.available))])  
      warning(paste(unique(wmod), collapse = "-")," in 'stkObs.models' is not an internal FLBEIA stock Observation model. If you want to use create.obs.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = immediate)
    }  
      
    if(!all(indObs.models %in% indObs.models.available)){ 
      wmod <- unique(indObs.models[which(!(indObs.models %in% indObs.models.available))])  
      warning(paste(unique(wmod), collapse = "-")," in 'indObs.models' is not an internal FLBEIA stock Observation model. If you want to use create.obs.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = immediate)
    }  
    
    # Generate the general structure
    k1 <- 1
    for(st in 1:length(stksnames)){
        res[[st]] <- vector('list', 2)
        names(res[[st]]) <- c("stkObs", "indObs")
        res[[st]][['stkObs']] <- list()
        res[[st]][['stkObs']][['stkObs.model']] <- stkObs.models[st]
        
        if(n.stks.inds[st] > 0){
            res[[st]][['indObs']] <- vector('list', n.stks.inds[st])
        
            for(id in 1:n.stks.inds[st]){
                res[[st]][['indObs']][['indObs.model']] <- indObs.models[k1]
                k1 <- k1 + 1 
            }
        }
    }
    
    # Add the function specific elements by stock/index.
    k1 <- 1  
    for(st in 1:length(stksnames)){  # fleet level functions; stkObs models
        
        stkObsmodcreator <- paste('create', stkObs.models[st],  'ctrl', sep = '.')

            
        res[[st]] <- eval(call(stkObsmodcreator, res = res[[st]], stkname = stksnames[st], largs = extra.args))
  
        if(n.stks.inds[st] > 0){
            for(id in 1:n.stks.inds[st]){ # index level functios:  indObs models.  
                
                indObsmodcreator <- paste('create', indObs.models[k1], 'ctrl', sep = '.')
            
                res[[st]][[id]] <- eval(call(indObsmodcreator, res = res[[st]][[id]], stkname = st, indname = id, largs = extra.args))
  
            
                k1 <- k1 + 1
            }
        }
    }
    
    return(res) 
} 

#-------------------------------------------------------------------------------
#                       ** create.perfectObs.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.perfectObs.ctrl <- function(resst,stkname,largs) return(resst)

#-------------------------------------------------------------------------------
#                       ** create.NoObsStock.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.NoObsStock.ctrl <- function(resst,stkname,largs) return(resst)

#-------------------------------------------------------------------------------
#                       ** create.NoObsIndex.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.NoObsIndex.ctrl <- function(resstid,stkname, indname, largs) return(resstid)

#-------------------------------------------------------------------------------
#                       ** create.bioInd.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.bioInd.ctrl <- function(resstid, stkname, indname, largs) return(resstid)

#-------------------------------------------------------------------------------
#                       ** create.ageInd.ctrl **
# error.ages - NULL or arra(na,na,ny, it)
#              o if is null => identity matrix, no error in aging.
#-------------------------------------------------------------------------------
create.ageInd.ctrl <- function(resstid,stkname, indname, largs){ 

    error.ages <- largs[[paste(error.ages,stkname, indname, sep = ".")]]
    flq.stk    <-  largs[[paste(flq,stkname, sep = ".")]]
    na <- dim(flq.stk)[1]
    ny <- dim(flq.stk)[2]
    it <- dim(flq.stk)[6]
    
    # if is null => identity matrix, no error in aging.
    if (is.null(error.ages)) {
        error.ages <- array(0, dim = c(na, na, ny, it), dimnames = list(dimnames(flq.stk)[[1]], 
            dimnames(flq.stk)[[1]], dimnames(flq.stk)[[2]], dimnames(flq.stk)[[6]]))
        for (a in 1:na) error.ages[a, a, , ] <- 1
    }
    
    resstid[['error.ages']] <- error.ages
    
return(resstid)
}

#-------------------------------------------------------------------------------
#                       ** create.age2ageDat.ctrl **
# flq.stkname MUST be provided to give 
# error.ages => default: identity matrix, no error in aging.
#
#-------------------------------------------------------------------------------
create.age2ageDat.ctrl <- function(resstid,stkname, indname, largs){ 

    flq.stk    <- largs[[paste('flq',stkname, sep = ".")]]
    
    
    if(is.null(flq.stk)) stop("You MUST provide 'flq.",stkname,"' object with '", stkname, "' stock's shape to be able to create 'age2ageDat' control object structure.")
    
    # collapse season dimension.
    flq.stk <- seasonSums(flq.stk)
    # extract main dimensions.
    na <- dim(flq.stk)[1]
    ny <- dim(flq.stk)[2]
    it <- dim(flq.stk)[6]
    
    # identity matrix, no error in aging.
    error.ages <- array(0, dim = c(na, na, ny, it), dimnames = list(dimnames(flq.stk)[[1]], 
                        dimnames(flq.stk)[[1]], dimnames(flq.stk)[[2]], dimnames(flq.stk)[[6]]))
    for (a in 1:na) error.ages[a, a, , ] <- 1
    
    # No error in any of the variables => FLQ = 1 for all.
    varia.mwgt <- varia.dwgt  <- varia.mort  <- varia.fec <- varia.ltot <- varia.dtot <- FLQuant(1, dimnames = dimnames(flq.stk))
    TAC.ovrsht <- FLQuant(1, dim = c(1, dim(flq.stk)[2],1,1,1,dim(flq.stk)[6]), dimnames = list(quant = 'all', year = dimnames(flq.stk)[[2]], iter = dimnames(flq.stk)[[6]]))
    
    resstid[['error.ages']] <- error.ages
    resstid[['varia.mwgt']] <- varia.mwgt
    resstid[['varia.dwgt']] <- varia.dwgt
    resstid[['varia.fec']]  <- varia.fec
    resstid[['varia.ltot']] <- varia.ltot
    resstid[['varia.dtot']] <- varia.dtot
    resstid[['TAC.ovrsht']] <- TAC.ovrsht
    
    return(resstid)
}


#-------------------------------------------------------------------------------
#                       ** create.age2agePop.ctrl **
#   Equal to create.age2ageDat.ctrl
#-------------------------------------------------------------------------------
create.age2agePop.ctrl <- function(resstid,stkname, indname, largs){ 

    flq.stk    <- largs[[paste(flq,stkname, sep = ".")]]
    na <- dim(flq.stk)[1]
    ny <- dim(flq.stk)[2]
    it <- dim(flq.stk)[6]
    
    # identity matrix, no error in aging.
    error.ages <- array(0, dim = c(na, na, ny, it), dimnames = list(dimnames(flq.stk)[[1]], 
                        dimnames(flq.stk)[[1]], dimnames(flq.stk)[[2]], dimnames(flq.stk)[[6]]))
    for (a in 1:na) error.ages[a, a, , ] <- 1
    
    # No error in any of the variables => FLQ = 1 for all.
    varia.ntot <- varia.mwgt <- varia.dwgt  <- varia.mort  <- varia.fec <- varia.ltot <- varia.dtot <- FLQuant(1, dimnames = dimnames(flq.stk))
    TAC.ovrsht <- FLQuant(1, dim = c(1, dim(flq.stk)[2],1,1,1,dim(flq.stk)[6]), dimnames = list(quant = 'all', year = dimnames(flq.stk)[[1]], iter = dimnames(flq.stk)[[6]]))
    
    resstid[['error.ages']] <- error.ages
    resstid[['varia.ntot']] <- varia.ntot
    resstid[['varia.mwgt']] <- varia.mwgt
    resstid[['varia.dwgt']] <- varia.dwgt
    resstid[['varia.fec']]  <- varia.fec
    resstid[['varia.ltot']] <- varia.ltot
    resstid[['varia.dtot']] <- varia.dtot
    resstid[['TAC.ovrsht']] <- TAC.ovrsht
    
    return(resstid)
}


#-------------------------------------------------------------------------------
#                       ** create.bio2bioDat.ctrl **
#-------------------------------------------------------------------------------
create.bio2bioDat.ctrl <- function(resstid,stkname, indname, largs){ 

    flq.stk    <- largs[[paste(flq,stkname, sep = ".")]]
    na <- dim(flq.stk)[1]
    ny <- dim(flq.stk)[2]
    it <- dim(flq.stk)[6]
    

    # No error in any of the variables => FLQ = 1 for all.
    varia.ltot <- varia.dtot <- TAC.ovrsht <- FLQuant(1, dimnames = dimnames(flq.stk)) 

    resstid[['varia.ltot']] <- varia.ltot
    resstid[['varia.dtot']] <- varia.dtot
    resstid[['TAC.ovrsht']] <- TAC.ovrsht
    
    return(resstid)
}


#-------------------------------------------------------------------------------
#                       ** create.bio2bioPop.ctrl **
#  Almost Equal to create.bio2bioDat.ctrl
#-------------------------------------------------------------------------------
create.bio2bioPop.ctrl <- function(resstid,stkname, indname, largs){ 
    
    flq.stk    <- largs[[paste(flq,stkname, sep = ".")]]
    na <- dim(flq.stk)[1]
    ny <- dim(flq.stk)[2]
    it <- dim(flq.stk)[6]
    

    # No error in any of the variables => FLQ = 1 for all.
    varia.btot <- varia.ltot <- varia.dtot <- TAC.ovrsht <- FLQuant(1, dimnames = dimnames(flq.stk))
       
    resstid[['varia.btot']] <- varia.btot
    resstid[['varia.ltot']] <- varia.ltot
    resstid[['varia.dtot']] <- varia.dtot
    resstid[['TAC.ovrsht']] <- TAC.ovrsht
    
    return(resstid)
}


#-------------------------------------------------------------------------------
#                       ** create.age2bioDat.ctrl **
#   Equal to create.bio2bioDat.ctrl
#-------------------------------------------------------------------------------
create.age2bioDat.ctrl <-  create.bio2bioDat.ctrl 


#-------------------------------------------------------------------------------
#                       ** create.age2bioPop.ctrl **
#   Equal to create.bio2bioPop.ctrl
#-------------------------------------------------------------------------------
create.age2bioPop.ctrl <- create.bio2bioPop.ctrl


