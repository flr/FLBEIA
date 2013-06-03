#-------------------------------------------------------------------------------
#                    ** create.observation.ctrl **
#
# Dorleta García - Azti Tecnalia
# 29/05/2013 11:34:54
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - **  ** : character vector with covariables names
# - **  ** : characted vector with the same length as cvrsnames with the process model followed by each of the covariables. 
#         the first element correspond with the process model of the first covariable in cvrsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCovar' is used for **all** the fleets.    

create.obs.ctrl <- function(stknames, n.stks.inds, stks.indsnames, stkObs.models = NULL, indObs.models = NULL,...){

    res        <- vector('list', length(stknames))
    names(res) <- stknames
    extra.args <- list(...)
    
    if(is.null(stkObs.models)) stkObs.models <- rep('perfectObs', length(stknames))
    if(is.null(indObs.models)) indObs.models <- rep('NoObsIndex', length(stks.indsnames))
    
    # Generate the general structure
    k1 <- 1
    for(st in 1:length(stknames)){
        res[[st]] <- vector('list', 2)
        names(res[[st]]) <- c("stkObs", "indObs")
        res[[st]][['stkObs']] <- list()
        res[[st]][['stkObs']][['stkObs.model']] <- stkObs.models[st]
        
        res[[st]][['indObs']] <- vector('list', n.stks.inds[st])
        
        for(id in 1:n.stks.inds[st]){
            res[[st]][['indObs']][['indObs.model']] <- indObs.models[k1]
            k1 <- k1 + 1 
        }
    }
    
    # Add the function specific elements by stock/index.
    k1 <- 1  
    for(st in 1:length(stknames)){  # fleet level functions; stkObs models
        
        stkObsmodcreator <- paste('create', stkObs.models[st],  'ctrl', sep = '.')

            
        res[[st]] <- call(stkObsmodcreator, res = res[[st]], stkname = st, largs = extra.args)
  
        for(id in 1:n.stks.inds[st]){ # index level functios:  indObs models.  
                
            indObsmodcreator <- paste('create', indObs.models[k1], 'ctrl', sep = '.')
            
            res[[st]][[id]] <- call(indObsmodcreator, res = res[[st]][[id]], stkname = st, indname = id, largs = extra.args)
  
            
            k1 <- k1 + 1
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
# error.ages => identity matrix, no error in aging.
#-------------------------------------------------------------------------------
create.age2ageDat.ctrl <- function(resstid,stkname, indname, largs){ 

    flq.stk    <- largs[[paste(flq,stkname, sep = ".")]]
    na <- dim(flq.stk)[1]
    ny <- dim(flq.stk)[2]
    it <- dim(flq.stk)[6]
    
    # identity matrix, no error in aging.
    error.ages <- array(0, dim = c(na, na, ny, it), dimnames = list(dimnames(flq.stk)[[1]], 
                        dimnames(flq.stk)[[1]], dimnames(flq.stk)[[2]], dimnames(flq.stk)[[6]]))
    for (a in 1:na) error.ages[a, a, , ] <- 1
    
    # No error in any of the variables => FLQ = 1 for all.
    varia.mwgt <- varia.dwgt  <- varia.mort  <- varia.fec <- varia.ltot <- varia.dtot <- FLQuant(1, dimnames = dimnames(flq.stk))
    TAC.ovrsht <- FLQuant(1, dim = c(1, dim(flq.stk)[2],1,1,1,dim(flq.stk)[6]), dimnames = list(quant = 'all', year = dimnames(flq.stk)[[1]], iter = dimnames(flq.stk)[[6]]))
    
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


