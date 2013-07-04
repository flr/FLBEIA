#-------------------------------------------------------------------------------
#                    ** create.covars.ctrl **
#
# Dorleta Garc?a - Azti Tecnalia
# 29/05/2013 14:04:15
#-------------------------------------------------------------------------------
#
#   :: ARGUMENTS ::
#
# - ** stksnames ** : character vector with stocks names
# - **  ** : characted vector with the same length as cvrsnames with the process model followed by each of the covariables. 
#         the first element correspond with the process model of the first covariable in cvrsnames, the second with the second and so on.
#         The default is NULL in which case 'fixedCovar' is used for **all** the fleets.    

create.advice.ctrl <- function(stksnames, HCR.models = NULL, HCR.ctrls = NULL,...){
    
    HCR.models.available <- c('fixedAdvice','annualTAC','IcesHCR','ghlHCR','annexIVHCR', 'froeseHCR')

    res        <- vector('list', length(stksnames))
    names(res) <- stksnames
    extra.args <- list(...)
    
    if(is.null(HCR.models)) HCR.models <- rep('fixedAdvice', length(stksnames))
    else{ 
        if(length(HCR.models) != length(stksnames)) stop("'HCR.models' must be NULL or must have the same length as stknames'")
        if(!all(HCR.models %in% HCR.models.available)){ 
            wmod <- unique(HCR.models[which(!(HCR.models %in% HCR.models.available))])  
            warning(paste(unique(wmod), collapse = "-")," in 'HCR.models' is not an internal FLBEIA covariables model. If you want to use create.covars.ctrl you must create, ", paste('create', paste(unique(wmod), collapse = ', ') ,'ctrl', sep = ".")," function.", immediate. = immediate)
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
#                        *** create.annexIVHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.fixedAdvice.ctrl <- function(resst,stkname, largs) return(resst)
    

#------------------------------------------------------------------------------#
#                        *** create.annualTAC.ctrl ***
#-------------------------------------------------------------------------------
create.annualTAC.ctrl <- function(resst,stkname, largs){

    resst <- c(resst, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = TRUE, 
                fwd.ctrl = NULL, sr = list(params = NULL, model = 'gmean', years = NULL),
                growth.years = NULL, advice = "TAC")
                
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
        
    resst$fwd.ctrl <- fwdControl(data.frame(year = c(0, 1, 1, 1),  val = c(1, Ftarget.stk, NA, NA), quantity = c( 'f', 'f', 'f', 'catch'),
                     min = c(NA,NA,0.9, 0.85), max  = c(NA,NA,1.1,1.15), rel.year = c(-1,NA,0, 0)))
    
   return(resst)
}

#------------------------------------------------------------------------------#
#                        *** create.IcesHCR.ctrl  ***
#-------------------------------------------------------------------------------
create.IcesHCR.ctrl <- function(resst,stkname, largs){

    resst <- c(resst, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = TRUE, 
                ref.pts = NULL, intermediate.year = 'Fsq', sr = list(params = NULL, model = 'gmean', years = NULL),
                growth.years = NULL, advice = "TAC")
                
    
    ref.pts.stk <- largs[[paste("ref.pts",stkname, sep = ".")]]
    
    
    
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
        
    if(!is.matrix(ref.pts.stk) | !all(c('Blim', 'Btrigger', 'Fmsy') %in% rownames(ref.pts.stk)))   stop(paste("ref.pts",stkname,sep = "."), " must be a matrix with dimension 3x(numb. of iterations) and rownames = c('Blim', 'Btrigger', 'Fmsy')")
    
    if(!is.null(largs$iter))  if(largs$iter != dim(ref.pts.stk)[2]) stop("Number of iterations in 'ref.pts.", stkname, "' must be equal to the iterations specified in 'iter' argument." )
 
    resst$ref.pts <- ref.pts.stk
    
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