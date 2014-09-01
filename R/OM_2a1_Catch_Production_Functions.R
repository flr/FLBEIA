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
# Dorleta GarcYYYa
# Created: 03/06/2011 10:35:14
# Changed: 03/06/2011 10:35:18
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#  CobbDouglasBio: E[it], B[it], q.m,efs.m,alpha.m,beta.m :: [mt,it] 
#       The functions _works_ with iterations
#-------------------------------------------------------------------------------
CobbDouglasBio   <- function(E,N, wl.m, wd.m, q.m,efs.m,alpha.m,beta.m,...)  # dga: aYYYado a como argumento.
                {
    Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m
    N   <- matrix(N,dim(efs.m)[2], byrow = T)
    C.m <-  q.m*(Ef*efs.m)^alpha.m*(N*(ret.m*wl.m + (1-ret.m)*wd.m))^beta.m  #
        
    C <-  colSums(C.m)

    return(catch =  C)
}


#-------------------------------------------------------------------------------
#  CobbDouglasBio.effort Cr[1], B[1], q.m,efs.m,alpha.m,beta.m :: [mt]       
#       The function does _not_work_ with iterations
#-------------------------------------------------------------------------------
CobbDouglasBio.effort   <- function(Cr,N, wl.m, wd.m,q.m,efs.m,alpha.m,beta.m,ret.m, restriction = 'catch',...){

    fObj <- function(E.f,Cr,N, wl.m, wd.m, q.m,efs.m,alpha.m,beta.m,ret.m, restriction){
        if(restriction == 'catch') C.m <- q.m*(E.f*efs.m)^alpha.m*(ret.m*(c(N)*wl.m)^beta.m + (1-ret.m)*(c(N)*wd.m)^beta.m)  # if restriction = catch (=> the restriction is catch not landings. )
        else C.m <- q.m*(E.f*efs.m)^alpha.m*(ret.m*(c(N)*wl.m))^beta.m
        return(Cr - sum(C.m))
    }
    
    # set upper limit
    X <- 10^(0:100)
    fobjX <- abs(sapply(X, fObj, Cr = Cr, N = N, wl.m = wl.m, wd.m = wd.m, q.m = q.m, efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, restriction = restriction, ret.m = ret.m))
    upl <- X[which(fobjX != Inf)[length(which(fobjX != Inf))]]
            
    NomEff <- uniroot(fObj,interval=c(0,upl),Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m,  q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, restriction = restriction, ret.m = ret.m)$root

    return(effort =  NomEff)
}
    


#-------------------------------------------------------------------------------
#  CobbDouglasAge :: E[it], B[na,nu,it], efs.m[mt,it], q.m,alpha.m,beta.m :: [mt,na,nu,it] 
#-------------------------------------------------------------------------------

CobbDouglasAge   <- function(E,N, wl.m, wd.m, ret.m,q.m,efs.m,alpha.m,beta.m,...){

    dimq <- dim(q.m)
    
    Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m      # [mt,it]
    Ef  <- array(Ef, dim = c(dim(Ef), dimq[2:3]))
    Ef  <- aperm(Ef, c(1,3:4,2))   # [mt,na,nu,it] 
   
    N <- array(N, dim = c(dim(N), dimq[1]))
    N <- aperm(N, c(4,1:3))      # [mt,na,nu,it]
    
    C.m <- q.m*(Ef*efs.m)^alpha.m*(N*(ret.m*wl.m + (1-ret.m)*wd.m))^beta.m # [mt,na,nu,it]
        
    C <-  apply(C.m, 4,sum)

    return(catch =  C)
}

#-------------------------------------------------------------------------------
#  CobbDouglasAge.Effort :: Cr[1], B[na,nu], efs.m[mt], q.m,alpha.m,beta.m :: [mt,na,nu] 
#-------------------------------------------------------------------------------

CobbDouglasAge.effort   <- function(Cr,N,wl.m, wd.m, ret.m, q.m,efs.m,alpha.m,beta.m, restriction = 'catch',...){
 
    dimq <- dim(q.m)  # [mt,na,nu,1]
    
    N <- array(N[drop = TRUE], dim(q.m)[2:4]) # [na,nu,1]
        
    N <- array(N, dim = c(dim(N), dimq[1]))  # [na,1,nu,1,1,1,mt]
    N <- aperm(N, c(4,1:3))                   # [mt,na,nu,it]
    
    efs.m <- array(efs.m, dim = dimq)

    fObj <- function(E.f,Cr,N,wd.m, wl.m, q.m,efs.m,alpha.m,beta.m, ret.m,restriction){
        # if catch = TRUE (=> the restriction is catch not landings. )
        if(restriction == 'catch') Ca.m <- q.m*(E.f*efs.m)^alpha.m*(N*(ret.m*wl.m+ (1-ret.m)*wd.m))^beta.m
        else  Ca.m <- q.m*(E.f*efs.m)^alpha.m*(N*ret.m*wl.m)^beta.m
        
                return(Cr - sum(Ca.m))
    }

    NomEff <- uniroot(fObj,interval=c(0,1e100),Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, restriction = restriction, ret.m = ret.m)$root

    return(effort =  NomEff)
}