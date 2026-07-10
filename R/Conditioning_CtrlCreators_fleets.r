#-------------------------------------------------------------------------------
#                    ** create.fleets.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 28/05/2013 10:42:07
#-------------------------------------------------------------------------------
#' fleets.ctrl object creator
#' 
#' It creates the fleets.ctrl object to be used in the call to the main function FLBEIA.
#' 
#'  
#   :: ARGUMENTS ::
#
#' @param fls character vector with fleet names
#' @param fls.mets.stksnames character vector containing the names of the stocks caught by each fleet–metier combination. The names attribute identifies the corresponding fleet, metier, and stock, using the format fleet::metier::stock.
#' @param catch.threshold if(NULL) => 0.9 for all the stocks (NULL is the default)  
#'                   else it must be an FLQuant with dim = c(nstks,ny,1,ns,nit)
#' @param seasonal.share an FLQuant with dimension [num. fleets, num. years, 1,num. seasons, 1, num. iterations] with elements between 0 and 1 to indicate how the quota of each fleet is distributed along seasons. The sum along seasons (seasonSums) must return an FLQuant with all elements equal to 1.  
#' @param effort.models  characted vector with the same length as fls with the effort model followed by each of the fleet. 
#'         the first element correspond with the effort model of the first fleet in fls, the second with the second and so on.
#'         The default is NULL in which case 'fixedEffort' is used for **all** the fleets.    
#' @param capital.models  characted vector with the same length as fls with the capital model followed by each of the fleet. 
#'         the first element correspond with the capital model of the first fleet in fls, the second with the second and so on.
#'         The default is NULL in which case 'fixedCapital' is used for **all** the fleets. 
#' @param catch.models  characted vector with the same length as sum(n.fls.stks) with the catch model followed by each of the fleet for each stock. 
#'         the first element correspond with the catch model of the first fleet in fls and the first stock in fls.stksnames, the second with the second and so on.
#'         The default is NULL in which case 'CobbDouglasAge' is used for **all** the fleets. 
#' @param price.models  characted vector with the same length as sum(n.fls.stks) with the price model followed by each of the fleet for each stock. 
#'         the first element correspond with the price model of the first fleet in fls and the first stock in fls.stksnames, the second with the second and so on.
#'         The default is NULL in which case 'fixedPrice' is used for **all** the fleets. 
#' @param flq An FLQuant to give structure to the FLQuants to be used within the function, 
#'         the dimension and dimnames in 'year', 'season' and 'iter' will be used to create the necessary FLQuants. 
#' @param ... Any extra arguments necessary in the model specific creators. '...' are extracted using 'list(...)', this generates a named list with the extra arguments.
#'        To assure the correct functioning the extra arguments must have a name, for example, elas = FLQuant(1,dimnames = DimsNms).
# 
# THIS TEXT IS ONLY RELEVANT FOR PEOPLE THAT GENERATE ITS OWN FUNCTIONS AND CREATORS TO RUN fleets.om
#               
#         and the creators for specific models must have the arguments 'res' and largs: 
#         \itemize{ 
#            \item res the element with the general structure created in the basice creator of fleets.ctrl object 
#            \item larg A named list with the elements necessary within the creator. For example if we use elas = FLQuant(1,dimnames = DimsNms) 
#                 in the call to the main function within the specific creator this element will be used as 'largs$elas'.
#           }         
#' @return A list of lists with the basic structure of the fleets.ctrl object.
#' 
 
create.fleets.ctrl <- function(fls, fls.mets.stksnames, catch.threshold = NULL,  seasonal.share = NULL, 
                                effort.models = NULL, capital.models = NULL, catch.models = NULL, price.models = NULL, flq, ...){
 
    effort.models.available  <- c('fixedEffort', 'SMFB', 'SSFB', 'MaxProfit', 'MaxProfitSeq')
    catch.models.available   <- c('CobbDouglasAge', 'CobbDouglasBio', 'Baranov')
    price.models.available   <- c('fixedPrice', 'elasticPrice', 'elasticPriceAge')
    capital.models.available <- c('fixedCapital', 'SCD')
    
    ## create hierarquical information fleet-metier-stock
    ids <- names(fls.mets.stksnames)
    
    fleet  <- sub("::.*", "", ids)
    metier <- sub("^[^:]+::(.*)::[^:]+$", "\\1", ids)
    stock  <- sub(".*::", "", ids)
    
    # n.fls.mets.stks
    fleet.metier <- paste(fleet, metier, sep = "::")
    fleet.metier.order <- unique(fleet.metier)
    n.fls.mets.stks <- table(
      factor(fleet.metier, levels = fleet.metier.order)
    )
    n.fls.mets.stks <- as.vector(n.fls.mets.stks)
    names(n.fls.mets.stks) <- fleet.metier.order
    
    # fls.metsnames
    idx <- !duplicated(fleet.metier)
    fls.metsnames <- setNames(
      metier[idx],
      fleet.metier[idx]
    )
    
    # fls.stksnames
    fleet.stock <- paste(fleet, stock, sep = "::")
    idx <- !duplicated(fleet.stock)
    fls.stksnames <- setNames(
      stock[idx],
      fleet.stock[idx]
    )
    
    # n.fls.stks
    idx <- !duplicated(fleet.stock)
    n.fls.stks <- table(
      factor(fleet[idx], levels = unique(fleet))
    )
    n.fls.stks <- as.vector(n.fls.stks)
    names(n.fls.stks) <- unique(fleet)
    
    
    # n.flts.mets
    idx <- !duplicated(fleet.metier)
    n.fls.mets <- table(
      factor(fleet[idx], levels = unique(fleet))
    )
    n.fls.mets <- as.vector(n.fls.mets)
    names(n.fls.mets) <- unique(fleet)
    
    nfls <- length(fls) 
    nmets <- length(fls.metsnames)
    res  <- vector('list', nfls + 2)
    names(res) <- c('catch.threshold', 'seasonal.share', fls)
    stknms <- unique(fls.stksnames)
    
    extra.args <- list(...) # a named list with the extra arguments in the call to the function, it'll contain the arguments for specific  creators.
  
    # Give default values to NULL elements.
    if(is.null(catch.threshold)){ res[['catch.threshold']] <- FLQuant(0.9, dimnames = list(stock = stknms, year = dimnames(flq)[['year']], 
                                                                        season = dimnames(flq)[['season']], iter = dimnames(flq)[['iter']]))
    }else{
        res[['catch.threshold']] <- catch.threshold
    }
    
    # Seasonal share.
    if(is.null(seasonal.share)){
        res[['seasonal.share']]        <- vector('list', length(unique(fls.stksnames)))
        names(res[['seasonal.share']]) <- unique(fls.stksnames)
        Dims <- c(list(fleet = fls), dimnames(flq)[2:6])
        for(i in 1: length(res[['seasonal.share']])) res[['seasonal.share']][[i]] <- FLQuant(1/dim(flq)[4], dimnames = Dims)
    }
    else{
        res[['seasonal.share']] <- seasonal.share
    }
    
    if(is.null(effort.models))  effort.models   <- rep('fixedEffort', nfls)                                                                      
    if(is.null(capital.models)) capital.models  <- rep('fixedCapital', nfls)  
    if(is.null(catch.models))   catch.models    <- rep('CobbDouglasAge', sum(n.fls.stks)) 
    if(is.null(price.models))   price.models    <- rep('fixedPrice', length(fls.mets.stksnames))  
    
    # check that all flq-s differ only in first (quant) dimension.
    test.flqs <- lapply(c('flq', names(extra.args)[grep(pattern = 'flq', names(extra.args))]), 
                        function(x) get(x))
    if(length(test.flqs)>1) if(!equal.flq.Dimnames(test.flqs, 2:5)) stop("All the input 'FLquant's must share 'year', 'unit', 'season', 'area' and 'iter' dimensions.")

    
    # Check foo.models
    # effort models
    if(length(effort.models) != nfls) stop("'effort.models' must be NULL or must have the same length as fls'")
    if(!all(effort.models %in% effort.models.available)){ 
        wmod <- effort.models[which(!(effort.models %in% effort.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'effort.models' is not an internal FLBEIA effort model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    # catch models.
    if(length(catch.models) != sum(n.fls.stks)) stop("'catch.models' must be NULL or must have the same length as stknames'")
    if(!all(catch.models %in% catch.models.available)){ 
        wmod <- catch.models[which(!(catch.models %in% catch.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'catch.models' is not an internal FLBEIA catch model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    # price models.
    if(length(price.models) != length(fls.mets.stksnames)) stop("'price.models' must be NULL or must have the same length as fleet-metier-stknames combination'")
    if(!all(price.models %in% price.models.available)){ 
        wmod <- price.models[which(!(price.models %in% price.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'price.models' is not an internal FLBEIA price model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    # capital models.
    if(length(capital.models) != nfls) stop("'capital.models' must be NULL or must have the same length as fls'")
    if(!all(capital.models %in% capital.models.available)){ 
        wmod <- capital.models[which(!(capital.models %in% capital.models.available))]  
        warning(paste(unique(wmod), collapse = ', ')," in 'capital.models' is not an internal FLBEIA capital model. If you want to use create.fleets.ctrl you must create, ", paste(paste('create', unique(wmod) ,'ctrl', sep = "."), collapse = ", ")," function.")
    }
    
    
    # Some check to test if declared number of stocks coincide with that used in the name declaration-
    # if(length(nflsts.mtrs) != nfls) stop("The length of 'fls' and 'nfls.mtrs' must coincide")
    # if(sum(nfls.mtrs) != length(fls.mtrsnames)) stop("The total number of metiers declared in 'nfls.mtrs' must be equal to the length of 'fls.mtrsnames'.")
    if(sum(n.fls.stks) != length(fls.stksnames)) stop("The total number of stocks declared in 'n.fls.stks' must be equal to the length of 'fls.stksnames'.")
    if(length(effort.models) != nfls) stop("The length of 'effort.models' must coincide with the number of fleets")
    if(length(capital.models) != nfls) stop("The length of 'capital.models' must coincide with the number of fleets")
    
    names(effort.models)  <- fls
    names(capital.models) <- fls
    names(n.fls.stks)    <- fls
    names(price.models)   <- fls.mets.stksnames
    names(catch.models)   <- fls.stksnames
    
    # Create the fleet/stock and/or fleet/metier/stock GENERAL structure in res. 
    k1 <- 1
    k2 <- 1
    k3 <- 1
    k4 <- 1
    k5 <- 1
    for(f in fls){  # fleet structure
        
        res[[f]]        <- vector('list', n.fls.stks[f] + 2) # minimum length
        names(res[[f]]) <- c('effort.model', 'capital.model', fls.stksnames[k1:(k1 + n.fls.stks[f] - 1)])
        
        res[[f]][['effort.model']]  <- unname(effort.models[f])
        res[[f]][['capital.model']] <- unname(capital.models[f])
        
        
        for(st in fls.stksnames[k1:(k1+n.fls.stks[f]-1)]){ # stock structure.  
            # res[[f]][[st]] <- vector('list', 1) # minimun length
            # names(res[[f]][[st]])  <- 'catch.model'
            # res[[f]][[st]][['catch.model']] <- unname(catch.models[k2])
            res[[f]][[st]] <- list(catch.model = catch.models[[k2]])
            k2 <- k2 + 1
        }
        k1 <- k1 + n.fls.stks[f]
        
        for(met in fls.metsnames[k3:(k3+n.fls.mets[f]-1)]){ # metier structure.  

          for(st in fls.mets.stksnames[k4:(k4+n.fls.mets.stks[paste0(f, "::", met)]-1)]){ # stock structure.  
            res[[f]][[met]][[st]] <- list(price.model = price.models[[k5]])
            k5 <- k5 + 1
          }
        k4 <- k4 + n.fls.mets.stks[paste0(f, "::", met)]
        }
        k3 <- k3 + n.fls.mets[f]
    }
    
    # Add the function specific elements by fleet/stock and/or fleet/metier/stock.
    k1 <- 1
    k2 <- 1
    k3 <- 1
    
    for(f in fls){  # fleet level functions; effort and capital models
      
        
        effmodcreator <- paste('create', effort.models[f],  'ctrl', sep = '.')
        capmodcreator <- paste('create', capital.models[f], 'ctrl', sep = '.')
            
        res[[f]] <- eval(call(effmodcreator, resf = res[[f]], fltname = f, largs = extra.args))
        res[[f]] <- eval(call(capmodcreator, resf = res[[f]], fltname = f, largs = extra.args))
                
        for(st in fls.stksnames[k1:(k1+n.fls.stks[f]-1)]){ # stock level functios: catch models.  
                
            catchmodcreator <- paste('create', catch.models[k1], 'ctrl', sep = '.')
            
            res[[f]][[st]] <- eval(call(catchmodcreator, resfst = res[[f]][[st]], fltname = f, stkname = st, largs = extra.args))
            
            k1 <- k1 + 1
        }
        for(met in fls.metsnames[k2:(k2+n.fls.mets[f]-1)]){ # metier-stock level functios: price models.  
          
          for(st in fls.mets.stksnames[k3:(k3+n.fls.mets.stks[paste0(f, "::", met)]-1)]){ # metier-stock level functios: price models.  
          pricemodcreator <- paste('create', price.models[k3], 'ctrl', sep = '.')
          
          res[[f]][[met]][[st]] <- eval(call(pricemodcreator, resfmetst = res[[f]][[met]][[st]], fltname = f, stkname = st, metname = met, largs = extra.args))
          
          k3 <- k3 + 1
          }
          k2 <- k2 + 1
        }
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
    
    rest      <- ifelse(is.null(largs[[paste('restriction', fltname,sep='.')]]), 'catch', largs[[paste('restriction', fltname,sep='.')]])
    effrest   <- ifelse(is.null(largs[[paste('effort.restr', fltname,sep='.')]]), 'max', largs[[paste('effort.restr', fltname,sep='.')]])

    
    if(is.na(effrest)) warning("Effort restriction in (effort.rest argument) is missing for fleet, ", fltname, ", NA used, you MUST fill it otherwise SMFB will not work.")
    
    resf[['restriction']]  <- rest 
    resf[['effort.restr']] <- effrest
    

    
    return(resf)
}

#-------------------------------------------------------------------------------
#                       ** create.MaxProfit.stkCnst.ctrl **
# extra args:  restriction.fltname, effort.restr.fltname
#-------------------------------------------------------------------------------
create.MaxProfit.ctrl <- function(resf,fltname,largs){
  
  # if NULL set default values, the 3rd stock in resf list, it can happen that is not an stock name.
  effort.restr   <- ifelse(is.null(largs[[paste('effort.restr', fltname,sep='.')]]),  NA, largs[[paste('effort.restr', fltname,sep='.')]])
  
  rest      <- ifelse(is.null(largs[[paste('restriction', fltname,sep='.')]]), NA, largs[[paste('restriction', fltname,sep='.')]])
  
  effort.range <- largs[[paste("effort.range", fltname, sep = ".")]]
  
  if(is.na( effort.restr )) warning("Effort constraint in (effort.restr argument) is missing for fleet, ",  fltname, ", NA used, you MUST fill it otherwise MaxProfit will not work.")
  if(is.na( rest )) warning("Catch/landing restriction in (srestriction argument) is missing for fleet, ",  fltname, ", NA used, you MUST fill it otherwise MaxProfit will not work.")
  
  resf[['effort.restr']] <- effort.restr
  resf[['restriction']]  <- rest 
  resf[["effort.range"]] <- effort.range
  
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
create.fixedPrice.ctrl <- function(resfmetst, fltname, metname, stkname, largs) return(resfmetst)


#-------------------------------------------------------------------------------
#                   ** create.elasticPrice.ctrl **
# extra args: flq.stkname
#-------------------------------------------------------------------------------
create.elasticPrice.ctrl <- function(resfmetst, fltname, metname, stkname, largs){
    
    flq.stk <- largs[[paste('flq', stkname,sep = ".")]]
    resfmetst[["pd.elsa"]]  <- flq.stk
    resfmetst[["pd.La0"]]   <- flq.stk
    resfmetst[["pd.Pa0"]]   <- flq.stk
    resfmetst[["pd.total"]] <- TRUE
    resfmetst[["type"]]     <- 1 # by default
    
    warning( "You have selected 'elasticPrice' price model for fleet, '", fltname, "',' metier", metname, "', and stock, '", stkname, 
             "'.\n Thus, you have to fill 'pd.elsa', 'pd.La0' and 'pd.Pa0' FLQuants in fleets.ctrl[['", fltname,"']][['", metname,"']][['", stkname,"']].\n",
             "'By default, the 'elasticPrice' model uses type 1. Select a different type if you want to apply another price model.\n ", sep = "")

    return(resfmetst)
}

#-------------------------------------------------------------------------------
#                   ** create.CobbDouglasAge.ctrl **
#                   ** create.CobbDouglasBio.ctrl **
#                   ** create.Baranov.ctrl **
#-------------------------------------------------------------------------------
create.CobbDouglasAge.ctrl <- create.CobbDouglasBio.ctrl <- create.Baranov.ctrl <- function(resfst, fltname, stkname, largs) return(resfst)

