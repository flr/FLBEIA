#-------------------------------------------------------------------------------
#           Auxiliary functions to ease coding of the complex functions
#
#       - unit2age()
#       - age2unit()
#       - unit2age() and age2unit() one inverse of the other one.
#       - dataCobbDoug()
# 
# Dorleta Garcia - 05/08/2010 12:12:19
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# unit2age(FLQuant[na,ny,nu,ns,1,it]) => array[na*nu,ny,ns,it]
#-------------------------------------------------------------------------------

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
setMethod("iter", signature(object="NULL"),
	  function(object, iter)
	  {
	         
		return(object)
	  }
) # }}}



 




