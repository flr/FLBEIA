#-------------------------------------------------------------------------------
#                    ** create.covars.ctrl **
#
# Dorleta García - Azti Tecnalia
# 29/05/2013 10:43:04
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - ** cvrsnames ** : character vector with covariables names
# - ** cvrs.models ** : characted vector with the same length as cvrsnames with the process model followed by each of the covariables. 
#         the first element correspond with the process model of the first covariable in cvrsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCovar' is used for **all** the fleets.    

create.covars.ctrl <- function(cvrsnames, cvrs.models = NULL, flq,...){

    res        <- vector('list', length(cvrsnames))
    names(res) <- cvrsnames
    extra.args <- list(...)
    
    if(is.null(cvrs.models)) cvrs.models <- rep('fixedCovar', length(cvrsnames))
    
    # Generate the general structure
    for(cv in 1:length(cvrsnames)){
        res[[cv]] <- list()
        res[[cv]][['process.model']] <- cvrs.models[cv] 
    }
    
    # Add process.model specific arguments.
    
    for(cv in 1:length(cvrsnames)){
        
        processmodcreator <- paste('create', process.models[cv],  'ctrl', sep = '.')
        res[[cv]] <- call(processmodcreator, res = res[[cv]], cvrname = cv, args = extra.args)
    }
    
    return(res) 
} 

#-------------------------------------------------------------------------------
#                       ** create.fixedCovar.ctrl **
#-------------------------------------------------------------------------------
create.fixedCovar.ctrl <- function(rescv,cvrname,largs) return(rescv)
