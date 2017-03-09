#-------------------------------------------------------------------------------
#                    ** create.assess.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 29/05/2013 14:04:15
#-------------------------------------------------------------------------------
#
#' assess.ctrl object creator
#' 
#' It creates the assess.ctrl object to be used in the call to the main function FLBEIA.
#' 
#
#   :: ARGUMENTS ::
#
#' @param stksnames A vector with the name of the stocks in the OM.
#' @param assess.models A character vector of the same length as stksnames with the name of the model used to obtaine the perceived population in the MP.
#' @param assess.ctrls A list of the same length as stksnames with the arguments needed to fit the assessment model.
# @param immediate logical, indicating if the warnings should be output immediately.
#' @param ... any extra arguments necessary in the model specific creators. '...' are extracted using 'list(...)', this generates a named list with the extra arguments.
#'        To assure the correct functioning the extra arguments must have a name.
#' 
#' @return A list of lists with the basic structure of the assess.ctrl object.
#-------------------------------------------------------------------------------


create.assess.ctrl <- function(stksnames, assess.models = NULL, assess.ctrls = NULL,...){

    assess.models.available <- c('NoAssessment', 'FLXSA') 
  
    res        <- vector('list', length(stksnames))
    names(res) <- stksnames
    extra.args <- list(...)
    
    if(is.null(assess.models)) assess.models <- rep('NoAssessment', length(stksnames))
    else{ 
      if(length(assess.models) != length(stksnames)) stop("'assess.models' must be NULL or must have the same length as stknames'")
      if(!all(assess.models %in% assess.models.available)){ 
        wmod <- unique(assess.models[which(!(assess.models %in% assess.models.available))])  
        warning(paste(unique(wmod), collapse = "-")," in 'assess.models' is not an internal FLBEIA covariables model. If you want to use create.covars.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = TRUE)
      }}
    
    
    # Generate the general structure
    for(st in 1:length(stksnames)){
        res[[st]] <- list()
        res[[st]][['assess.model']] <- assess.models[st] 
    }
    
    # Add assessment controls, they must be in the call to the function, otherwise they are not created here.
    
    for(st in 1:length(stksnames)){
        
        assessmodcreator <- paste('create', assess.models[st],  'ctrl', sep = '.')
        res[[st]] <- eval(call(assessmodcreator, res = res[[st]], stkname = stksnames[st], largs = extra.args))
    }
    
    return(res) 
} 


#-------------------------------------------------------------------------------
#                       ** create.NoAssessment.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.NoAssessment.ctrl <- function(resst,stkname,largs) return(resst)

#-------------------------------------------------------------------------------
#                       ** create.FLXSA.ctrl **
# No extra arguments needed
#-------------------------------------------------------------------------------
create.FLXSA.ctrl <- function(resst,stkname,largs){ 
    resst <- c(resst, control = FLXSA::FLXSA.control())
    warning("Default values have been used to create 'FLXSA.control' for stock '", stkname, "' change it by hand if the values are not appropiate.")
    
    return(resst)}


