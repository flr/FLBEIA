#-------------------------------------------------------------------------------
#                   CATCH PRODUCTION FUNCTIONS 
#  'catch.model' argument in fleets.ctrl[[flnm]] object
#  Catch production function defines the relationship between effort and catch.
#
#   - CobbDouglasBio - Catch from Cobb-Douglass production function aggregated in biomass
#   - CobbDouglasAge - Catch from Cobb-Douglass production function by age.
#
#   - CobbDouglasBio.effort - Effort from Cobb-Douglass production function aggregated in biomass
#   - CobbDouglasAge.effort - Effort from Cobb-Douglass production function by age.
#
# Dorleta García
# Created: 03/06/2011 10:35:14
# Changed: 03/06/2011 10:35:18
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#  CobbDouglasBio: E[it], B[it], q.m,efs.m,alpha.m,beta.m :: [mt,it] 
#       The functions _works_ with iterations
#-------------------------------------------------------------------------------
CobbDouglasBio   <- function(E,B,q.m,efs.m,alpha.m,beta.m,...)  # dga: añado a como argumento.
                {
    Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m
    B   <- matrix(B,dim(efs.m)[2], byrow = T)
    C.m <- q.m*Ef^alpha.m*B^beta.m
        
    C <-  colSums(C.m)

    return(catch =  C)
}


#-------------------------------------------------------------------------------
#  CobbDouglasBio.effort Cr[1], B[1], q.m,efs.m,alpha.m,beta.m :: [mt]       
#       The functions does _not_work_ with iterations
#-------------------------------------------------------------------------------
CobbDouglasBio.effort   <- function(Cr,B,q.m,efs.m,alpha.m,beta.m,...){

    fObj <- function(E.f,Cr,B, q.m,efs.m,alpha.m,beta.m){
        C.m <- q.m*(E.f*efs.m)^alpha.m*B^beta.m
        return(Cr - sum(C.m))
    }
    
    # set upper limit
    X <- 10^(0:100)
    fobjX <- abs(sapply(X, fObj, Cr = Cr, B = B, q.m = q.m, efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m))
    upl <- X[which(fobjX != Inf)[length(which(fobjX != Inf))]]
            
    NomEff <- uniroot(fObj,interval=c(0,upl),Cr=Cr,B=B, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m)$root

    return(effort =  NomEff)
}
    


#-------------------------------------------------------------------------------
#  CobbDouglasAge :: E[it], B[na,nu,it], efs.m[mt,it], q.m,alpha.m,beta.m :: [mt,na,nu,it] 
#-------------------------------------------------------------------------------

CobbDouglasAge   <- function(E,Ba,q.m,efs.m,alpha.m,beta.m,...){

    dimq <- dim(q.m)
    
    Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m      # [mt,it]
    Ef  <- array(Ef, dim = c(dim(Ef), dimq[2:3]))
    Ef  <- aperm(Ef, c(1,3:4,2))   # [mt,na,nu,it] 
   
    Ba <- array(Ba, dim = c(dim(Ba), dimq[1]))
    Ba <- aperm(Ba, c(4,1:3))      # [mt,na,nu,it]
    
    C.m <- q.m*(Ef^alpha.m)*(Ba^beta.m)   # [mt,na,nu,it]
        
    C <-  apply(C.m, 4,sum)

    return(catch =  C)
}

#-------------------------------------------------------------------------------
#  CobbDouglasAge.Effort :: Cr[1], B[na,nu], efs.m[mt], q.m,alpha.m,beta.m :: [mt,na,nu] 
#-------------------------------------------------------------------------------

CobbDouglasAge.effort   <- function(Cr,Ba,q.m,efs.m,alpha.m,beta.m,...){
 
    dimq <- dim(q.m)  # [mt,na,nu,1]
    
    Ba <- array(Ba[drop = TRUE], dim(q.m)[2:4]) # [na,nu,1] 
        
    Ba <- array(Ba, dim = c(dim(Ba), dimq[1]))    # [na,1,nu,1,1,1,mt]
    Ba <- aperm(Ba, c(4,1:3))
    
    efs.m <- array(efs.m, dim = dimq)

    fObj <- function(E.f,Cr,Ba, q.m,efs.m,alpha.m,beta.m){
        Ca <- q.m*(E.f*efs.m)^alpha.m*(Ba^beta.m)
     #   cat('Ca_1', Ca, '\n')
  #      Ca <- ifelse(Ca>Ba, Ba*0.95, Ca) 
     #   cat('Ca_2', Ca, '\n')
        C.m <- sum(Ca)
        
        return(Cr - sum(C.m))
    }

    NomEff <- uniroot(fObj,interval=c(0,1e100),Cr=Cr,Ba=Ba, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m)$root

    return(effort =  NomEff)
}