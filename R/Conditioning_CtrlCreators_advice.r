#-------------------------------------------------------------------------------
#                    ** create.advice.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 29/05/2013 14:04:15
# Modified
# 9/11/2016  Agurtzane Urtizberea
#-------------------------------------------------------------------------------
#
#' advice.ctrl object creator
#' 
#' It creates the advice.ctrl object to be used in the call to the main function FLBEIA.
#' 
#
#   :: ARGUMENTS ::
#
#' @param stksnames A vector with the name of the stocks in the OM.
#' @param HCR.models A character vector of the same length as stksnames with the name of the HCR used to generate the management advice.
#' @param first.year The first year in which advice is calculated.
#' @param last.year The last year in which advice is calculated.
# @param immediate logical, indicating if the warnings should be output immediately.
#' @param ... any extra arguments necessary in the HCR specific creators. '...' are extracted using 'list(...)', this generates a named list with the extra arguments.
#'        To assure the correct functioning the extra arguments must have a name.
#' 
#' @return A list of lists with the basic structure of the advice.ctrl object.
#-------------------------------------------------------------------------------

create.advice.ctrl <- function(stksnames, HCR.models = NULL, ...){
    
    HCR.models.available <- c('fixedAdvice','annualTAC','IcesHCR','ghlHCR','annexIVHCR', 'froeseHCR', 'F2CatchHCR', 'neaMAC_ltmp', 'aneHCR_JD')

    res        <- vector('list', length(stksnames))
    names(res) <- stksnames
    extra.args <- list(...)
    
    if(is.null(HCR.models)) HCR.models <- rep('fixedAdvice', length(stksnames))
    else{ 
        if(length(HCR.models) != length(stksnames)) stop("'HCR.models' must be NULL or must have the same length as stknames'")
        if(!all(HCR.models %in% HCR.models.available)){ 
            wmod <- unique(HCR.models[which(!(HCR.models %in% HCR.models.available))])  
            warning(paste(unique(wmod), collapse = "-")," in 'HCR.models' is not an internal FLBEIA covariables model. If you want to use create.covars.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = TRUE)
        }}
    
    
    # Generate the general structure
    for(st in 1:length(stksnames)){
        res[[st]] <- list()
        res[[st]][['HCR.model']] <- HCR.models[st] 
    }
    
    # Add advicement controls, they must be in the call to the function, otherwise they are not created here.
    
    for(st in 1:length(stksnames)){
        
        HCRmodcreator <- paste('create', HCR.models[st],  'ctrl', sep = '.')
        res[[st]] <- eval(call(HCRmodcreator, res = res[[st]], stkname = stksnames[st], largs = extra.args))
    }
    
    return(res) 
} 


#------------------------------------------------------------------------------#
#                        *** create.fixedAdvice.ctrl  ***
#-------------------------------------------------------------------------------
create.fixedAdvice.ctrl <- function(resst,stkname, largs) return(resst)
    

#------------------------------------------------------------------------------#
#                        *** create.annualTAC.ctrl ***
#-------------------------------------------------------------------------------
create.annualTAC.ctrl <- function(resst,stkname, largs){
   
    first.yr <- largs$first.yr  
    last.yr <- largs$last.yr

    resst <- c(resst, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = TRUE, 
                fwd.ctrl = NULL,AdvCatch=NULL,  
                growth.years = NULL,first.yr=NULL,last.yr=NULL)
                
    resst$sr <- list(params = NULL, model = 'geomean', years = NULL)

    AdvCatch.stk <- largs[[paste("AdvCatch", stkname, sep = ".")]]
    
    cat("--------------------- NOTE ON ADVICE ------------------------------------------------------------------------------\n")            
    cat("A default control for 'annualTAC' HCR has been created for stock, ", stkname,".\n")
    cat("Fadvice = Ftarget subject to: \n * 10% maximum variation in F and,\n * 15% maximum variation in TAC.\n")
    cat("In the intermediate year fishing mortality equal to F statu quo.\n")
    cat("For recruitment or population growth in biomass a geometric mean of historic time series estimates will be used.\n")
    cat("Average of last 3 years used for biological parameters and fishing mortality.\n")
    
    Ftarget.stk <- largs[[paste("Ftarget",stkname, sep = ".")]]
    if(is.null(Ftarget.stk)){ 
        warning("Ftarget for stock, '", stkname,"' has not been specified in argument: ", paste("Ftarget",stkname,sep = "."), ".\n - A fwd.ctrl element with empty Ftarget will be created. FILL IT BY HAND!!!!")
        Ftarget.stk <- NA
    }
    cat("--------------------------------------------------------------------------------------------------------------------\n") 
    
    if (is.null(first.yr)| is.null(last.yr)) {
      stop("first.yr (first year with historic data) and last.yr (last year of projection) must be defined")
      cat("------------------------------------------------------------------------------\n")
    } 
    
    if (is.null(AdvCatch.stk)) {
      AdvCatch.stk <- rep(FALSE,length(first.yr:last.yr))
      names(AdvCatch.stk) <- c(first.yr:last.yr)
      warning("Advice of ", stkname, " is FALSE by default, so the advice is given in terms of landings 
              ", immediate. = TRUE)
      cat("------------------------------------------------------------------------------\n")
    }    
    
    resst$fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1, 1, 1),  val = c(1, Ftarget.stk, NA, NA), quantity = c( 'f', 'f', 'f', 'catch'),
                     min = c(NA,NA,0.9, 0.85), max  = c(NA,NA,1.1,1.15), rel.year = c(-1,NA,0, 0)))

    resst$AdvCatch <- AdvCatch.stk
    
   return(resst)
}

#------------------------------------------------------------------------------#
#                        *** create.IcesHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.IcesHCR.ctrl <- function(resst,stkname, largs){

    first.yr <- largs$first.yr  
    last.yr <- largs$last.yr
    
    resst <- c(resst, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = TRUE, 
                ref.pts = NULL, AdvCatch=NULL, intermediate.year = 'Fsq',
                growth.years = NULL, first.yr=NULL,last.yr=NULL)
    resst$sr <- list(params = NULL, model = 'geomean', years = NULL)
                
    
    ref.pts.stk <- largs[[paste("ref.pts",stkname, sep = ".")]]
    AdvCatch.stk <- largs[[paste("AdvCatch", stkname, sep = ".")]]
    
    
    cat("--------------------- NOTE ON ADVICE ------------------------------------------------------------------------------\n")            
    cat("A default control for 'IcesHCR' HCR has been created for stock, ", stkname,".\n")
    cat("In the intermediate year fishing mortality equal to F statu quo.\n")
    cat("For recruitment or population growth in biomass a geometric mean of historic time series estimates will be used.\n")
    cat("Average of last 3 years used for biological parameters and fishing mortality.\n")
    
    
    if(is.null(ref.pts.stk)){
        it <- ifelse(is.null(largs$iter), 1, largs$iter)
        warning("Reference points for stock, '", stkname,"' have not been specified in argument: ", paste("ref.pts",stkname,sep = "."), ". \n -  A ref.pts element with empty reference points has been created. FILL IT BY HAND!!!!", immediate. = TRUE)
        if(is.null(it))  warning("iter argument is missing, iter = 1 will be used in the creation of ref.pts element, correct it if necessary.")
        ref.pts.stk <- matrix(NA, 3,it, dimnames = list( c('Blim', 'Btrigger', 'Fmsy'), 1:it))
        cat("------------------------------------------------------------------------------\n") 
    }
    
    if (is.null(first.yr)| is.null(last.yr)) {
      stop("first.yr (first year with historic data) and last.yr (last year of projection) must be defined")
      cat("------------------------------------------------------------------------------\n")
    } 
    
    if (is.null(AdvCatch.stk)) {
      AdvCatch.stk <- rep(FALSE,length(first.yr:last.yr))
      names(AdvCatch.stk) <- c(first.yr:last.yr)
      warning("Advice of ", stkname, " is FALSE by default, so the advice is given in terms of landings 
              ", immediate. = TRUE)
      cat("------------------------------------------------------------------------------\n")
    }       
        
    if(!is.matrix(ref.pts.stk) | !all(c('Blim', 'Btrigger', 'Fmsy') %in% rownames(ref.pts.stk)))   stop(paste("ref.pts",stkname,sep = "."), " must be a matrix with dimension 3x(numb. of iterations) and rownames = c('Blim', 'Btrigger', 'Fmsy')")
    
    if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
 
    resst$ref.pts <- ref.pts.stk
    resst$AdvCatch <- AdvCatch.stk

    return(resst)
}



#------------------------------------------------------------------------------#
#                        *** create.annexIVHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.annexIVHCR.ctrl <- function(resst,stkname, largs){

    resst <- c(resst, index = 1, ref.pts = NULL, type = 2)
                
    cat("--------------------- NOTE ON ADVICE --------------------------------\n")            
    cat("A default control for 'annexIVHCR' HCR has been created for stock, ", stkname,"\n")
    cat("The first index will be used to apply the HCR using type = 2. \n")

    
    ref.pts.stk <- largs[[paste("ref.pts",stkname, sep = ".")]]
    
    if(is.null(ref.pts.stk)){
        it <- ifelse(is.null(largs$iter), 1, largs$iter)
        warning("Reference points for stock, '", stkname,"' have not been specified in argument: ", paste("ref.pts",stkname,sep = "."), ". \n -  A ref.pts element with alpha = 0.20 and beta = 0.15 reference points will be created.\n - If you don't agree with the values FILL THEM BY HAND!!!!")
        if(is.null(it))  warning("iter argument is missing, iter = 1 will be used in the creation of ref.pts element, correct it if necessary.")
        ref.pts.stk <- matrix(NA, 2,it, dimnames = list( c('alpha', 'beta'), 1:it))
        ref.pts.stk[1,] <- 0.2
        ref.pts.stk[2,] <- 0.15
    }
      
    if(!is.matrix(ref.pts.stk) | !all(c('alpha', 'beta') %in% rownames(ref.pts.stk)))  stop(paste("ref.pts",stkname,sep = "."), " must be a matrix with dimension 2x(numb. of iterations) and rownames = c('alpha', 'beta')")
    
    cat("------------------------------------------------------------------------\n")  
    
    if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
    
    
    resst$ref.pts <- ref.pts.stk

   return(resst)
}



#------------------------------------------------------------------------------#
#                        *** create.ghlHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.ghlHCR.ctrl <- function(resst,stkname, largs){

    resst <- c(resst, ref.pts = NULL)
                
    cat("--------------------- NOTE ON ADVICE --------------------------------\n")            
    cat("A default control for 'ghlHCR' HCR has been created for stock, ", stkname,"\n")

    
    ref.pts.stk <- largs[[paste("ref.pts",stkname, sep = ".")]]
    
    if(is.null(ref.pts.stk)){
        it <- ifelse(is.null(largs$iter), 1, largs$iter)
        warning("Reference points for stock, '", stkname,"' have not been specified in argument: ", paste("ref.pts",stkname,sep = "."), ". \n -  A ref.pts element with alpha_0 = 2, alpha_1 = 1 and beta = 0.05 reference points will be created.\n - If you don't agree with the values FILL THEM BY HAND!!!!")
        if(is.null(it))  warning("iter argument is missing, iter = 1 will be used in the creation of ref.pts element, correct it if necessary.")
        ref.pts.stk <- matrix(NA, 3,it, dimnames = list( c('alpha_0', 'alpha_1', 'beta'), 1:it))
        ref.pts.stk[1,] <- 2
        ref.pts.stk[2,] <- 1
        ref.pts.stk[3,] <- 0.05
    }
    
    if(!is.matrix(ref.pts.stk) | !all(c('alpha_0', 'alpha_1', 'beta') %in% rownames(ref.pts.stk)))  stop(paste("ref.pts",stkname,sep = "."), " must be a matrix with dimension 3x(numb. of iterations) and rownames = c('alpha_0', 'alpha_1', 'beta')")

    cat("------------------------------------------------------------------------\n") 
    
    if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
    
         
    resst$ref.pts <- ref.pts.stk
    
   return(resst)
}


#------------------------------------------------------------------------------#
#                        *** create.froeseHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.froeseHCR.ctrl <- function(resst,stkname, largs){

    resst <- c(resst, ref.pts = NULL)
                
    cat("--------------------- NOTE ON ADVICE --------------------------------\n")            
    cat("A default control for 'froeseHCR' HCR has been created for stock, ", stkname,"\n")

    
    ref.pts.stk <- largs[[paste("ref.pts",stkname, sep = ".")]]
    
    if(is.null(ref.pts.stk)){
        it <- ifelse(is.null(largs$iter), 1, largs$iter)
        warning("Reference points for stock, '", stkname,"' have not been specified in argument: ", paste("ref.pts",stkname,sep = "."), ". \n -  A ref.pts element with alpha_0 = 1, alpha_1 = 5, beta = 0.91 and empty Bmsy and MSY reference points will be created.\n - FILL IN Bmsy and MSY BY HAND and If you don't agree with the rest of the values FILL THEM ALSO  BY HAND!!!!")
        if(is.null(it))  warning("iter argument is missing, iter = 1 will be used in the creation of ref.pts element, correct it if necessary.")
        ref.pts.stk <- matrix(NA, 5,it, dimnames = list( c('alpha_0', 'alpha_1', 'beta', 'Bmsy', 'MSY'), 1:it))
        ref.pts.stk[1,] <- 1
        ref.pts.stk[2,] <- 0.5
        ref.pts.stk[3,] <- 0.91
    }
    
    if(!is.matrix(ref.pts.stk) | !all(c('alpha_0', 'alpha_1', 'beta','Bmsy', 'MSY') %in% rownames(ref.pts.stk)))  stop(paste("ref.pts",stkname,sep = "."), " must be a matrix with dimension 5x(numb. of iterations) and rownames = c('alpha_0', 'alpha_1', 'beta', 'MSY' and 'Bmsy')")

    cat("------------------------------------------------------------------------\n") 
    
    if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
    
    resst$ref.pts <- ref.pts.stk
    
   return(resst)
}


#------------------------------------------------------------------------------#
#                        *** create.F2CatchHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.F2CatchHCR.ctrl <- function(resst,stkname, largs){
  
   first.yr <- largs$first.yr  
   last.yr <- largs$last.yr
  

    resst <- c(resst, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = TRUE, 
                ref.pts = NULL, AdvCatch=NULL, intermediate.year = 'Fsq',
                growth.years = NULL,first.yr=NULL,last.yr=NULL)
    resst$sr <- list(params = NULL, model = 'geomean', years = NULL)
                
    
    ref.pts.stk <- largs[[paste("ref.pts",stkname, sep = ".")]]
    AdvCatch.stk <- largs[[paste("AdvCatch", stkname, sep = ".")]]
    
    
    cat("--------------------- NOTE ON ADVICE ------------------------------------------------------------------------------\n")            
    cat("A default control for 'F2CatchHCR' HCR has been created for stock, ", stkname,".\n")
    cat("In the intermediate year fishing mortality equal to F statu quo.\n")
    cat("For recruitment or population growth in biomass a geometric mean of historic time series estimates will be used.\n")
    cat("Average of last 3 years used for biological parameters and fishing mortality.\n")
    
    
    if(is.null(ref.pts.stk)){
        it <- ifelse(is.null(largs$iter), 1, largs$iter)
        warning("Reference points for stock, '", stkname,"' have not been specified in argument: ", paste("ref.pts",stkname,sep = "."), ". \n -  A ref.pts element with empty reference points has been created. FILL IT BY HAND!!!!", immediate. = TRUE)
        if(is.null(it))  warning("iter argument is missing, iter = 1 will be used in the creation of ref.pts element, correct it if necessary.")
        ref.pts.stk <- matrix(NA, 1,it, dimnames = list( c('Ftarget'), 1:it))
        cat("------------------------------------------------------------------------------\n") 
    }
    if (is.null(first.yr)| is.null(last.yr)) {
      stop("first.yr (first year with historic data) and last.yr (last year of projection) must be defined")
      cat("------------------------------------------------------------------------------\n")
    } 
    if (is.null(AdvCatch.stk)) {
      AdvCatch.stk <- rep(FALSE,length(first.yr:last.yr))
      names(AdvCatch.stk) <- c(first.yr:last.yr)
      warning("Advice of ", stkname, " is FALSE by default, so the advice is given in terms of landings 
              ", 
              immediate. = TRUE)
      cat("------------------------------------------------------------------------------\n")
    }    
        
    if(!is.array(ref.pts.stk) | !all(c('Ftarget') %in% dimnames(ref.pts.stk)[[1]]))   stop(paste("ref.pts",stkname,sep = "."), " must be a array with dimension 1x(numb. projections years)x(numb. of iterations) and dimnames[[1]] = c('Ftarget')")
    
    if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
  
    resst$ref.pts <- ref.pts.stk
    resst$AdvCatch <- AdvCatch.stk
   return(resst)
}


#------------------------------------------------------------------------------#
#                        *** create.neaMAC_ltmp.ctrl  ***
#-------------------------------------------------------------------------------
create.neaMAC_ltmp.ctrl <- create.IcesHCR.ctrl


#------------------------------------------------------------------------------#
#                        *** create.aneHCRs.ctrl  ***
#-------------------------------------------------------------------------------
create.aneHCR_JD.ctrl <- function(resst,stkname, largs){
  
  first.yr <- largs$first.yr  
  last.yr <- largs$last.yr
  
  resst <- c(resst, ref.pts = NULL, AdvCatch=NULL, first.yr=NULL, last.yr=NULL) 
  
  ref.pts.stk         <- largs[[paste("ref.pts",stkname, sep = ".")]]
  TACs1.perc.stk      <- largs[[paste("TACs1.perc",stkname, sep = ".")]]
  tsurv.stk           <- largs[[paste("tsurv",stkname, sep = ".")]]
  cbbm.params.flq.ANE <- largs[[paste("cbbm.params.flq",stkname, sep = ".")]]
  
  AdvCatch.stk <- largs[[paste("AdvCatch", stkname, sep = ".")]]
  
  if(is.null(ref.pts.stk)){
    it <- ifelse(is.null(largs$iter), 1, largs$iter)
    warning("Reference points for stock, '", stkname,"' have not been specified in argument: ", paste("ref.pts",stkname,sep = "."), ". \n 
            -  A ref.pts element with empty reference points has been created. FILL IT BY HAND!!!!", immediate. = TRUE)
    if(is.null(it))  warning("iter argument is missing, iter = 1 will be used in the creation of ref.pts element, correct it if necessary.")
    ref.pts.stk <- matrix(NA, 7,it, dimnames = list( c('alpha','gamma','TACmin', 'TACmax', 'Btrig1', 'Btrig2', 'Btrig3'), 1:it))
    cat("------------------------------------------------------------------------------\n") 
  }
  
  if(!is.matrix(ref.pts.stk) | !all(c('alpha','gamma','TACmin', 'TACmax', 'Btrig1', 'Btrig2', 'Btrig3') %in% rownames(ref.pts.stk)))   
    stop(paste("ref.pts",stkname,sep = "."), " must be a matrix with dimension 7x(numb. of iterations) and rownames = c('alpha','gamma','TACmin', 'TACmax', 'Btrig1', 'Btrig2', 'Btrig3')")
  
  if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
  
  if (is.null(TACs1.perc.stk)) stop("Percentage of TAC captured in 1st season required for ", stkname, "', to be specified in 'TACs1.perc.", stkname,"' argument." )
  if (is.null(tsurv.stk)) stop("Moment of the survey required for ", stkname, "', to be specified in 'tsurv.", stkname,"' argument." )
  
  if (is.null(AdvCatch.stk)) {
    AdvCatch.stk <- rep(FALSE,length(first.yr:last.yr))
    names(AdvCatch.stk) <- c(first.yr:last.yr)
    warning("Advice of ", stkname, " is FALSE by default, so the advice is given in terms of landings 
            ", immediate. = TRUE)
    cat("------------------------------------------------------------------------------\n")
  } 
  
  cat("--------------------- NOTE ON ADVICE ------------------------------------------------------------------------------\n")            
  cat("Remember fill 'cbbm.params' for stock, ", stkname,"\n")
  
  resst$ref.pts     <- ref.pts.stk
  resst$TACs1.perc  <- TACs1.perc.stk
  resst$tsurv       <- tsurv.stk
  resst$cbbm.params <- list( G=cbbm.params.flq.ANE, M=cbbm.params.flq.ANE, S=cbbm.params.flq.ANE)
  
  index <- largs[[paste("index",stkname, sep = ".")]]
  if( is.null(index)) stop("Index required for ", stkname, "', to be specified in 'index.", stkname,"' argument." )
  resst$index <- index
  
  resst$AdvCatch <- AdvCatch.stk
  
  return(resst)
}

