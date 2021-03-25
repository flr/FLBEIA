#-------------------------------------------------------------------------------
#           Auxiliary functions to ease coding of the complex functions
#
#       - updateBiols()
#       - unit2age()
#       - age2unit()
#       - unit2age() and age2unit() one inverse of the other one.
#       - dataCobbDoug()
# 
# Dorleta Garcia - 05/08/2010 12:12:19
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# updateFLBiols(biols): Update the FLBiols object to the FLCore 2.6 
#-------------------------------------------------------------------------------

#' Update the FLBiols object to the FLCore 2.6
#' 
#' Updates an old FLBiols (where slots fec and mat from each FLBiol are FLQuants), into the new version of the class,
#' where slots fec and mat are of class \code{predictModel}.
#'
#' @param biols A FLBiols object.
#' 
#' @return A new FLBEIA output object with all the iterations joined. 
#' 
#' @seealso \code{\link{predictModel}}
#' 

updateFLBiols <- function(biols){
  
  res        <- vector('list', length(biols))
  names(res) <- names(biols)
  
  lapply(biols, function(x){
    mat0   <- x@fec
    mat0[] <- 1
    n0 <- x@n
    res   <- FLBiol(name = x@name,
                    desc = x@desc,
                    range = x@range,
                    n = x@n,
                    m = x@m,
                    wt = x@wt,
                    rec = predictModel(n0 = n0, model = ~n0[1,]),
                    mat = predictModel(mat = mat0, model = ~ mat),
                    fec = predictModel(fec = x@fec, model = ~ fec),
                    spwn = x@spwn) 
    res@range[4:5] <- c(dims(res@n)$minyear, dims(res@n)$maxyear)
    return(res)
  })
}

#-------------------------------------------------------------------------------
# unit2age(FLQuant[na,ny,nu,ns,1,it]) => array[na*nu,ny,ns,it]
#-------------------------------------------------------------------------------

#' Translates unit dimension in an FQuant into age dimension
#' 
#' This function transforms a FQuant with several 'unit' dimension into 
#' an array an unique 'unit' dimension and no 'area' dimension. Moving 'unit' to 'age' dimension.
#' This is usefull when we use the 'unit' dimension to store the different seasonal cohorts.
#'
#' @param Q A FLQuant.
#'               The object must be the output of FLBEIA function.
#'               
#' @return A 4-dimensional array with length c(d[1]*d[3], d[2], d[4], d[6]), 
#'         where d is the dimension of the FLQuant. 
#' 
#' @note The files must contain a single object (named as \code{objnam} value).
#' 

unit2age <- function(Q){

    Dim <- dim(Q)
    
    if(dim(Q)[3] == 1){
        A <- array(unname(Q[drop= TRUE]), dim = Dim[c(1,2,4,6)])
        return(A)
    }
        
    A   <- array(dim = c(Dim[1]*Dim[3], Dim[2], Dim[4], Dim[6]))
    
    k <- (0:(Dim[1]-1))*Dim[3] + 1

    for(ss in 1:Dim[4]){
        for(uu in 1:Dim[3]){
  #          cat('season: ',ss, ', unit:', uu,'\n')
            k1 <- k + (ss - uu)
            k2 <- ifelse(k1[1] <=0,2,1)
            A[k1[k1>0],,ss,] <- Q[k2:Dim[1],,uu,ss,,drop = TRUE]
        }
    }
    
   # convert to 0 the NAs in the plusgroups. (see next explanation)
   na <- ((Dim[1]*Dim[3]) - (Dim[4]-1))
   B <- A[-(1:na),,,,drop = FALSE] 
   B[is.na(B)]  <- 0
   A[-(1:na),,,]  <- B
   return(A)
         
}


#-------------------------------------------------------------------------------
# age2unit(array[na,ny,ns,it]) => FLQuant[na,ny,nu,ns,1,it])
#-------------------------------------------------------------------------------
# A: The array to transform
# Q: An FLQuant with the dimension and  dimnames we want for the resulting FLQ.

age2unit <- function(A,Q){
    Dim <- dim(Q) 
    B   <- array(0,dim = Dim)
    na  <- Dim[1]
    nu  <- Dim[3]   
    ns  <- Dim[4]
    
    if(dim(A)[1] == Dim[1]){   
        Q[] <- A    
        return(Q)    
    }
   
   # => nu = ns
    
    k <- (0:(na-1))*nu
    
    for(ss in 1:ns){
        for(uu in 1:nu){
   #          cat('season: ',ss, ', unit:', uu,'\n')
            # Which position has each unit in the age matrix?
            if(ss == nu) pos <- ss:1
            else pos <- c(ss:1,nu:(ss+1))
            
            pos <- which(pos == uu)
            k1 <- k + pos
            k0 <- 1:na
            
            if(uu>ss){
                k1 <- k1[-na]
                k0 <- 2:na
            }
            
            B[k0, , uu, ss, , ] <- A[k1, , ss,,drop=F ]
        } 
    }
    B <- FLQuant(B, dimnames = dimnames(Q))
    
    return(B)
}



## iter {{{
setMethod("iter", signature(obj="NULL"),
	  function(obj, iter)
	  {
	         
		return(obj)
	  }
) # }}}



#-------------------------------------------------------------------------------
# extractBRP(advice.ctrl): Extracts reference points from advice.ctrl object
#                          for calculating summary statistics with bioSum
#-------------------------------------------------------------------------------

#' Extracts reference points (specifically Bpa, Blim, Bmsy, Fpa, Flim and Fmsy) from advice.ctrl object,
#' If not available, then values are set to NA.
#' 
#' The output can be used for the brp argument of \code{bioSum} function.
#'
#' @param advice.ctrl A list with advice controls as the FLBEIA function argument advice.ctrl.
#' @param stkn Names of the stocks.
#' @param Btarget Named vector with the name of the target biological reference point for each element in stkn. 
#'                Default "Bmsy".
#' @param Ftarget Named vector with the name of the target fishing mortality reference point for each element in stkn. 
#'                Default "Fmsy".
# @inheritParams FLBEIA
#' 
#' @return A data frame with columns stock, iter and one colum per reference point with the value 
#'         of the biological reference points per stock and iteration. The used reference points are 
#'         Bpa, Blim, Bmsy, Fpa, Flim and Fmsy.
#' 
#' @seealso \code{\link{bioSum}}
#' 

#' @examples
#'\dontrun{
#'
#' library(FLBEIA)
#'
#' data(one)
#' extractBRP(oneAdvC, stkn = names(oneBio))
#' 
#' data(oneIt)
#' extractBRP(oneItAdvC, stkn = names(oneItBio))
#' 
#' data(multi)
#' extractBRP(multiAdvC, stkn = names(multiBio))
#'                           
#' # setting targets different to Bmsy and Fmsy
#' 
#' extractBRP(oneAdvC, stkn = names(oneBio),
#'            Btarget=setNames("Btrigger","stk1"), Ftarget=setNames("Ftrigger","stk1"))
#' 
#' Btarget <- setNames(c("Btrigger","Bmsy"),names(oneBio))
#' Ftarget <- setNames(c("Ftrigger","Fmsy"),names(oneBio))
#' 
#' extractBRP(multiAdvC, stkn = names(multiBio))
#' 
#' }

extractBRP <- function(advice.ctrl, stkn, Btarget=NULL, Ftarget=NULL) {
  
  
  refpts <- c("Bpa", "Blim", "Btarget", "Fpa", "Flim", "Ftarget")
  
  if (is.null(Btarget))
    Btarget <- setNames( rep("Bmsy", length(stkn)), stkn)
  
  if (is.null(Ftarget))
    Ftarget <- setNames( rep("Fmsy", length(stkn)), stkn)
  
  its <- lapply(advice.ctrl, function(y) dimnames(y$ref.pts)[[2]] %>% as.numeric()) %>% unlist() %>% unique()
  
  # check class
  
  if (!is.list(advice.ctrl)) stop("advice.ctrl argument must be a list")
  if (!is.character(stkn))  stop("stkn argument must be a character vector")
  if (!is.character(Btarget) | !all(names(Btarget)==stkn))  
    stop("Btarget must be a named character vector with one element per stkn")
  if (!is.character(Ftarget) | !all(names(Ftarget)==stkn))  
    stop("Ftarget must be a named character vector with one element per stkn")
  
  # check if there is one element per stkn
  
  if (any(!stkn %in% names(advice.ctrl))) 
    stop( paste0("Missing stocks in the list: ", paste(stkn[which(!stkn %in% names(advice.ctrl))], collapse = ", ")))
  
  # out <- setNames(data.frame(matrix(ncol = 8, nrow = 0)),c("stock", "iter", refpts))
  
  out <- data.frame(stock = character(), iter  = numeric(),
                    Bpa   = numeric(), Blim  = numeric(), Btarget  = numeric(), 
                    Fpa   = numeric(), Flim  = numeric(), Ftarget  = numeric())
  
  for (st in stkn) {
    
    refpts.st <- c("Bpa", "Blim", Btarget[[st]], "Fpa", "Flim", Ftarget[[st]])
    
    rp <- advice.ctrl[[st]]$ref.pts 
    
    if(is.null(rp) | !any(dimnames(rp)[[1]] %in% refpts.st)) {
      
      rp.df <- data.frame(stock = st, 
                          iter = its)
      
    } else {
      
      # check dimensions
      if (any(as.numeric(dimnames(rp)[[2]]) != its)) 
        stop( paste("Check iterations in reference points for",st))
      
      # extract values
      rp.df <- rp %>% t() %>% 
        as.data.frame() %>% 
        mutate(stock = st, 
               iter = as.numeric(dimnames(rp)[[2]])) %>% 
        select(stock, iter, any_of(refpts.st))
      
      rownames(rp.df) <- NULL
      
      if (Btarget[[st]] %in% names(rp.df))
        rp.df <- rp.df %>% dplyr::rename(Btarget = Btarget[[st]])
      
      if (Ftarget[[st]] %in% names(rp.df))
        rp.df <- rp.df %>% dplyr::rename(Ftarget = Ftarget[[st]])
      
    }
    
    out <- out %>% 
      bind_rows(rp.df)
    
  }
  
  
  return(out)
  
}


