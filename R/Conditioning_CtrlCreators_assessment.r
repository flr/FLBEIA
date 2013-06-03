#-------------------------------------------------------------------------------
#                    ** create.covars.ctrl **
#
# Dorleta García - Azti Tecnalia
# 29/05/2013 14:04:15
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - ** stksnames ** : character vector with stocks names
# - **  ** : characted vector with the same length as cvrsnames with the process model followed by each of the covariables. 
#         the first element correspond with the process model of the first covariable in cvrsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCovar' is used for **all** the fleets.    

create.assess.ctrl <- function(stksnames, assess.models = NULL, assess.ctrls = NULL){

    res        <- vector('list', length(stksnames))
    names(res) <- stksnames
    extra.args <- list(...)
    
    if(is.null(assess.models)) assess.models <- rep('NoAssessment', length(stksnames))
    
    # Generate the general structure
    for(st in 1:length(stksnames)){
        res[[st]] <- list()
        res[[st]][['assess.model']] <- assess.models[st] 
    }
    
    # Add assessment controls, they must be in the call to the function, otherwise they are not created here.
    
    for(st in 1:length(stksnames)){
        
        assessmodcreator <- paste('create', assess.models[st],  'ctrl', sep = '.')
        res[[st]] <- call(assessmodcreator, res = res[[st]], stkname = st, args = extra.args)
    }
    
    return(res) 
} 

