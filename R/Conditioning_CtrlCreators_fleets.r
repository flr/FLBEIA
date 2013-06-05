#-------------------------------------------------------------------------------
#                    ** create.fleets.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 28/05/2013 10:42:07
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - ** fltsnames ** : character vector with fleet names
# - ** n.flts.stks ** : numeric vector with the same length as fltsnames with the  declaration the stocks caugth by each of the fleets.
# - ** flts.stksnames **: character vector with length = sum(n.flts.stks), with the names of the stocks caught by the fleet, 
#        o the  first n.flts.stks[1] elements correspond to the stocks caught by the first fleet in fltsnames
#        o the  following n.flts.stks[2] elements correspond to the stocks caught by the second fleet in fltsnames
#        o and so on.
# - ** catch.threshold **: if(NULL) => 0.9 for all the stocks (NULL is the default)  
#                   else it must be an FLQuant with dim = c(nstks,ny,1,ns,nit)
# - ** effort.models **:  characted vector with the same length as fltsnames with the effort model followed by each of the fleet. 
#         the first element correspond with the effort model of the first fleet in fltsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedEffort' is used for **all** the fleets.    
# - ** capital.models **:  characted vector with the same length as fltsnames with the capital model followed by each of the fleet. 
#         the first element correspond with the capital model of the first fleet in fltsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCapital' is used for **all** the fleets. 
# - ** catch.models **:  characted vector with the same length as sum(n.flts.stks) with the catch model followed by each of the fleet for each stock. 
#         the first element correspond with the catch model of the first fleet in fltsnames and the first stock in flts.stksnames, the second with the second and so on.
#         The default is NULL in which case 'CobbDouglasAge' is used for **all** the fleets. 
# - ** price.models **:  characted vector with the same length as sum(n.flts.stks) with the price model followed by each of the fleet for each stock. 
#         the first element correspond with the price model of the first fleet in fltsnames and the first stock in flts.stksnames, the second with the second and so on.
#         The default is NULL in which case 'fixedPrice' is used for **all** the fleets. 
# - ** flq **: An flquant to give structure to the FLQuants to be used within the function, 
#         the dimensioan and dimnames in 'year', 'season' and 'iter' will be used to create the necessary FLQuants. 
# - ...: Any extra arguments necessary in the model specific creators. '...' are extracted using 'list(...)', this generates a named list with the extra arguments.
#        To assure the correct functioning the extra arguments must have a name, for example, elas = FLQuant(1,dimnames = DimsNms) and 
#        the creators for specific models must have to arguments 'res' and largs: 
#           o 'res': the element with the general structure created in the main function 
#           o 'larg': A named list with the elements necessary within the creator. For example if we use elas = FLQuant(1,dimnames = DimsNms) 
#                 in the call to the main function within the specific creator this element will be used as 'largs$elas'.
#                    

 
create.fleets.ctrl <- function(fltsnames,  n.flts.stks, flts.stksnames, catch.threshold = NULL,
                                effort.models = NULL, capital.models = NULL, catch.models = NULL, price.models = NULL, flq, ...){
    
    effort.models.available  <- c('fixedEffort', 'SMFB', 'SSFB', 'MaxProfit.stkCnst')
    catch.models.available   <- c('CobbDouglasAge', 'CobbDouglasBio')
    price.models.available   <- c('fixedPrice', 'elasticPrice')
    capital.models.available <- c('fixedCapital', 'SCD')
    
    nfls <- length(fltsnames) 
    res  <- vector('list', nfls + 1)
    names(res) <- c('catch.threshold', fltsnames)
    stknms <- unique(flts.stksnames)
    
    extra.args <- list(...) # a named list with the extra arguments in the call to the function, it'll contain the arguments for specific  creators.
  
    # Give default values to NULL elements.
    if(is.null(catch.threshold)) res[['catch.threshold']] <- FLQuant(0.9, dimnames = list(stock = stknms, year = dimnames(flq)[['year']], 
                                                                        season = dimnames(flq)[['season']], iter = dimnames(flq)[['iter']]))
    if(is.null(effort.models))  effort.models   <- rep('fixedEffort', nfls)                                                                      
    if(is.null(capital.models)) capital.models  <- rep('fixedCapital', nfls)  
    if(is.null(catch.models))   catch.models    <- rep('CobbDouglasAge', sum(n.flts.stks)) 
    if(is.null(price.models))   price.models    <- rep('fixedPrice', sum(n.flts.stks))  
    
    
    # Check foo.models
    # effort models
    if(length(effort.models) < nfls) stop("'effort.models' must be NULL or must have the same length as stknames'")
    if(!all(effort.models %in% effort.models.available)){ 
        wmod <- effort.models[which(!(effort.models %in% effort.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'effort.models' is not an internal FLBEIA effort model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    # catch models.
    if(length(catch.models) < sum(n.flts.stks)) stop("'catch.models' must be NULL or must have the same length as stknames'")
    if(!all(catch.models %in% catch.models.available)){ 
        wmod <- catch.models[which(!(catch.models %in% catch.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'catch.models' is not an internal FLBEIA catch model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    # price models.
    if(length(price.models) < sum(n.flts.stks)) stop("'price.models' must be NULL or must have the same length as stknames'")
    if(!all(price.models %in% price.models.available)){ 
        wmod <- price.models[which(!(price.models %in% price.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'price.models' is not an internal FLBEIA price model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    # capital models.
    if(length(capital.models) < nfls) stop("'capital.models' must be NULL or must have the same length as stknames'")
    if(!all(capital.models %in% capital.models.available)){ 
        wmod <- capital.models[which(!(capital.models %in% capital.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'capital.models' is not an internal FLBEIA capital model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    
    
    # Some check to test if declared number of stocks coincide with that used in the name declaration-
#    if(length(nflsts.mtrs) != nfls) stop("The length of 'fltsnames' and 'nflts.mtrs' must coincide")
#    if(sum(nflts.mtrs) != length(flts.mtrsnames)) stop("The total number of metiers declared in 'nflts.mtrs' must be equal to the length of 'flts.mtrsnames'.")
    if(sum(n.flts.stks) != length(flts.stksnames)) stop("The total number of stocks declared in 'n.flts.stks' must be equal to the length of 'flts.stksnames'.")
    if(length(effort.models) != nfls) stop("The length of 'effort.models' must coincide with the number of fleets")
    if(length(capital.models) != nfls) stop("The length of 'capital.models' must coincide with the number of fleets")
    
  
    names(effort.models)  <- fltsnames
    names(capital.models) <- fltsnames
    names(n.flts.stks)    <- fltsnames
    names(price.models)   <- flts.stksnames
    names(catch.models)   <- flts.stksnames
    
    # Create the fleet/metier/stock GENERAL structure in res. 
    k1 <- 1
    k2 <- 1
    for(f in fltsnames){  # fleet structure
        
        res[[f]]        <- vector('list', n.flts.stks[f] + 2) # minimum length
        names(res[[f]]) <- c('effort.model', 'capital.model', flts.stksnames[k1:(k1 + n.flts.stks[f] - 1)])
        
        res[[f]][['effort.model']]  <- unname(effort.models[f])
        res[[f]][['capital.model']] <- unname(capital.models[f])
        
        
        for(st in flts.stksnames[k1:(k1+n.flts.stks[f]-1)]){ # stock structure.  
            res[[f]][[st]] <- vector('list', 2) # minimun length
            names(res[[f]][[st]])  <- c('catch.model', 'price.model')
            res[[f]][[st]][['catch.model']] <- unname(catch.models[k2])
            res[[f]][[st]][['price.model']] <- unname(price.models[k2])
            k2 <- k2 + 1
        }
        k1 <- k1 + n.flts.stks[f]
    }
    
    # Add the function specific elements by fleet/stock.
    k1 <- 1
    jf <- 1
    for(f in fltsnames){  # fleet level functions; effort and capital models
      
        
        effmodcreator <- paste('create', effort.models[f],  'ctrl', sep = '.')
        capmodcreator <- paste('create', capital.models[f], 'ctrl', sep = '.')
            
        res[[f]] <- eval(call(effmodcreator, resf = res[[f]], fltname = f, largs = extra.args))
        res[[f]] <- eval(call(capmodcreator, resf = res[[f]], fltname = f, largs = extra.args))
                
        for(st in flts.stksnames[k1:(k1+n.flts.stks[jf]-1)]){ # stock level functios: catch and price models.  
                
            catchmodcreator <- paste('create', catch.models[k1], 'ctrl', sep = '.')
            pricemodcreator <- paste('create', price.models[k1], 'ctrl', sep = '.')
            
            res[[f]][[st]] <- eval(call(catchmodcreator, resfst = res[[f]][[st]], fltname = f, stkname = st, largs = extra.args))
            res[[f]][[st]] <- eval(call(pricemodcreator, resfst = res[[f]][[st]], fltname = f, stkname = st, largs = extra.args))
            
            k1 <- k1 + 1
        }
        jf <- jf + 1
    }
    
    return(res) 

} 

#-------------------------------------------------------------------------------
#                       ** create.fixedEffort.ctrl **
#-------------------------------------------------------------------------------
create.fixedEffort.ctrl <- function(resf,fltname,largs) return(resf)

#-------------------------------------------------------------------------------
#                       ** create.SMFB.ctrl **
# extra args: restriction.fltname, effort.restr.fltname
#-------------------------------------------------------------------------------
create.SMFB.ctrl <- function(resf, fltname,largs){

    # if NULL set default values
    rest    <- ifelse(is.null(largs[[paste('restriction', fltname,sep='.')]]), 'landings', largs[[paste('restriction', fltname,sep='.')]])
    effrest <- ifelse(is.null(largs[[paste('effort.restr', fltname,sep='.')]]), NA, largs[[paste('effort.restr', fltname,sep='.')]])
    
    if(is.null(largs[[paste('restriction', fltname,sep='.')]])) warning("Restriction in (effort.rest argument) is missing for fleet, ", fltname, ", 'landings' used,  if you want 'catch' restriction you MUST change it manually.")
    if(is.na(effrest)) warning("Effort restriction in (effort.rest argument) is missing for fleet, ", fltname, ", NA used, you MUST fill it otherwise SMFB will not work.")
    
    resf[['restriction']]  <- rest 
    resf[['effort.restr']] <- effrest
    return(resf)
}

#-------------------------------------------------------------------------------
#                       ** create.MaxProfit.stkCnst.ctrl **
# extra args:  stk.cnst.fltname
#-------------------------------------------------------------------------------
create.MaxProfit.stkCnst.ctrl <- function(resf,fltname,largs){
    
    # if NULL set default values, the 3rd stock in resf list, it can happen that is not an stock name.
    stk.cnst    <- ifelse(is.null(largs[[paste('stk.cnst', fltname,sep='.')]]),  NA, largs[[paste('stk.cnst', fltname,sep='.')]])

    if(is.na( stk.cnst )) warning("Stock constraint in (stk.cnst argument) is missing for fleet, ",  fltname, ", NA used, you MUST fill it otherwise MaxProfit.stkCnst will not work.")
    
    resf[['stk.cnst']] <- stk.cnst
 
    return(resf)
}

#-------------------------------------------------------------------------------
#                       ** create.fixedCapital.ctrl **
#-------------------------------------------------------------------------------
create.fixedCapital.ctrl <- function(resf,fltname,largs) return(resf)

#-------------------------------------------------------------------------------
#                       ** create.SCD.ctrl **
# extra args: NONE
# This function does not contain extra arguments in fleets.ctrl, the 
# extra information is stored in covars because they are system variables.
#-------------------------------------------------------------------------------
create.SCD.ctrl <- function(resf,fltname,largs) return(resf)


#-------------------------------------------------------------------------------
#                   ** create.fixedPrice.ctrl **
#-------------------------------------------------------------------------------
create.fixedPrice.ctrl <- function(resfst, fltname, stkname, largs) return(resfst)


#-------------------------------------------------------------------------------
#                   ** create.elasticPrice.ctrl **
# extra args: flq.stkname
#-------------------------------------------------------------------------------
create.elasticPrice.ctrl <- function(resfst, fltname, stkname, largs){
    
    flq.stk <- largs[[paste('flq', stkname,sep = ".")]]
    resfst[["pd.els"]]   <- flq.stk
    resfst[["pd.La0"]]   <- flq.stk
    resfst[["pd.Pa0"]]   <- flq.stk
    resfst[["pd.total"]] <- TRUE
    
    warning( "You have selected 'elasticPrice' price model for fleet, '", fltname, "', and stock, '", stkname, 
        "'.\n Thus, you have to fill 'pd.els', 'pd.La0' and 'pd.Pa0' FLQuants in fleets.ctrl[['", fltname,"']][['", stkname,"']].\n", sep = "")

    return(resfst)
}

#-------------------------------------------------------------------------------
#                   ** create.CobbDouglasAge.ctrl **
#                   ** create.CobbDouglasBio.ctrl **
#-------------------------------------------------------------------------------
create.CobbDouglasAge.ctrl <- create.CobbDouglasBio.ctrl <- function(resfst, fltname, stkname, largs) return(resfst)

