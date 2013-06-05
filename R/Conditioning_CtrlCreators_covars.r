#-------------------------------------------------------------------------------
#                    ** create.covars.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 29/05/2013 10:43:04
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - ** cvrsnames ** : character vector with covariables names
# - ** process.models ** : characted vector with the same length as cvrsnames with the process model followed by each of the covariables. 
#         the first element correspond with the process model of the first covariable in cvrsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCovar' is used for **all** the fleets.    

create.covars.ctrl <- function(cvrsnames, process.models = NULL, flq, immediate = FALSE,...){

    process.models.available <- 'fixedCovar'
  
    res        <- vector('list', length(cvrsnames))
    names(res) <- cvrsnames
    extra.args <- list(...)
    
    if(is.null(process.models)) process.models <- rep('fixedCovar', length(cvrsnames))
    else{ 
      if(length(process.models) < length(cvrsnames)) stop("'process.models' must be NULL or must have the same length as stknames'")
      if(!all(process.models %in% process.models.available)){ 
        wmod <- unique(process.models[which(!(process.models %in% process.models.available))])  
        warning(paste(unique(wmod), collapse = "-")," in 'process.models' is not an internal FLBEIA covariables model. If you want to use create.covars.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = immediate)
      }}
    
    
    # Generate the general structure
    for(cv in 1:length(cvrsnames)){
        res[[cv]] <- list()
        res[[cv]][['process.model']] <- process.models[cv] 
    }
    
    # Add process.model specific arguments.
    
    for(cv in 1:length(cvrsnames)){
        
        processmodcreator <- paste('create', process.models[cv],  'ctrl', sep = '.')
        res[[cv]] <- eval(call(processmodcreator, res = res[[cv]], cvrname = cv, largs = extra.args))
    }
    
    return(res) 
} 

#-------------------------------------------------------------------------------
#                       ** create.fixedCovar.ctrl **
#-------------------------------------------------------------------------------
create.fixedCovar.ctrl <- function(rescv,cvrname,largs) return(rescv)
