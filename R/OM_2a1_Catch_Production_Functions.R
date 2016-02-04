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
#       The function _works_ with iterations
# ONLY FOR BIOMASS DYNAMIC POPUYLATUIONS.-
#-------------------------------------------------------------------------------
CobbDouglasBio   <- function(E,N, wl.m, wd.m, q.m,efs.m,alpha.m,beta.m, ret.m, rho = 1 ,...)  # dga: aYYYado a como argumento.
                {
    Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m
    N   <- matrix(N,dim(efs.m)[1], dim(efs.m)[2],byrow = T)
    alpha.m <- matrix(alpha.m,dim(efs.m)[1],dim(efs.m)[2])
    wl.m <- matrix(wl.m,dim(efs.m)[1],dim(efs.m)[2])
    wd.m <- matrix(wd.m,dim(efs.m)[1],dim(efs.m)[2])
    beta.m <- matrix(beta.m,dim(efs.m)[1],dim(efs.m)[2])
    ret.m <- matrix(ret.m,dim(efs.m)[1],dim(efs.m)[2])
    q.m <- matrix(q.m,dim(efs.m)[1],dim(efs.m)[2])

    C.m <-  q.m*(Ef)^alpha.m*(N*(ret.m*wl.m + (1-ret.m)*wd.m))^beta.m  #
        
  #  C <-  colSums(C.m)


    Clim <- sweep(N*(ret.m*wl.m + (1-ret.m)*wd.m), 2, rho, "*")
    
    C.m <- ifelse(C.m < Clim,  C.m, Clim)
  
  #  C.m <- ifelse(C.m < rho*N*(ret.m*wl.m + (1-ret.m)*wd.m),  C.m, rho*N*(ret.m*wl.m + (1-ret.m)*wd.m))

    return(catch =  C.m)        # [nmt,it]
}


#-------------------------------------------------------------------------------
#  CobbDouglasBio.effort Cr[1], B[1], q.m,efs.m,alpha.m,beta.m :: [mt]       
#       The function does _not_work_ with iterations
#-------------------------------------------------------------------------------
CobbDouglasBio.effort   <- function(Cr,N, wl.m, wd.m,q.m,efs.m,alpha.m,beta.m,ret.m, rho = 1, restriction = 'catch',...){


    fObj <- function(E.f,Cr,N, wl.m, wd.m, q.m,efs.m,alpha.m,beta.m,ret.m, rho, restriction){


         C.m <- CobbDouglasBio(E = E.f, N = N, wl.m = wl.m, wd.m = wd.m, q.m = q.m, efs.m = efs.m,
                                 alpha.m = alpha.m, beta.m = beta.m, ret.m = ret.m, rho = rho)

        if(restriction == 'catch') C.m <- C.m  # if restriction = catch (=> the restriction is catch not landings. )
        else C.m <- matrix(ret.m,dim(efs.m)[1],dim(efs.m)[2])*C.m     # if restriction = landings
        return(Cr - sum(C.m))
    }
    
    # set upper limit
    X <- 10^(0:100)
    fobjX <- abs(sapply(X, fObj, Cr = Cr, N = N, wl.m = wl.m, wd.m = wd.m, q.m = q.m, efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, restriction = restriction, ret.m = ret.m, rho = rho))
    upl <- X[which(fobjX != Inf)[length(which(fobjX != Inf))]]

    Cinf <- CobbDouglasBio(E = upl,N=N, wl.m = wl.m, wd.m = wd.m, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, ret.m = ret.m, rho = rho)
    if((Cr - sum(Cinf))> 0) # Even with infinity effort it is not possible to catch the quota => return 'almost' infinity effort.
        return(effort = upl)


    NomEff <- uniroot(fObj,interval=c(0,upl),Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m,  q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, restriction = restriction, ret.m = ret.m, rho = rho)$root

    return(effort =  NomEff)
}
    


#-------------------------------------------------------------------------------
#  CobbDouglasAge :: E[it], B[na,nu,it], efs.m[mt,it], q.m,alpha.m,beta.m :: [mt,na,nu,it] 
# Res => C.m[mt,na,nu,it]
# OUTPUT: Catch at age in weight.
#-------------------------------------------------------------------------------

CobbDouglasAge   <- function(E,N, wl.m, wd.m, ret.m,q.m,efs.m,alpha.m,beta.m,rho = 1,...){

    dimq <- dim(q.m)
    
    Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m      # [mt,it]
    Ef  <- array(Ef, dim = c(dim(Ef), dimq[2:3]))
    Ef  <- aperm(Ef, c(1,3:4,2))   # [mt,na,nu,it] 
   
    N <- array(N, dim = c(dim(N), dimq[1]))
    N <- aperm(N, c(4,1:3))      # [mt,na,nu,it]
    
    W <- (ret.m*wl.m + (1-ret.m)*wd.m)

    C.m <- q.m*(Ef)^alpha.m*(N*W)^beta.m # [mt,na,nu,it]
    
    Clim <- sweep(W*N,4,rho,"*")

    C.m <- ifelse(C.m < Clim, C.m, Clim) # The truncation  of CobDoug is applied at metier level.

#    C.m <- ifelse(C.m < rho*W*N, C.m,rho*W*N) # The truncation  of CobDoug is applied at metier level.
    
 #   C <-  apply(C.m, 4,sum)

    return(catch =  C.m)  # [mt,na,nu,it]
}

#-------------------------------------------------------------------------------
#  CobbDouglasAge.Effort :: Cr[1], B[na,nu], efs.m[mt], q.m,alpha.m,beta.m :: [mt,na,nu] 
#-------------------------------------------------------------------------------

CobbDouglasAge.effort   <- function(Cr,N,wl.m, wd.m, ret.m, q.m,efs.m,alpha.m,beta.m, rho = 1, restriction = 'catch',...){
 


    fObj <- function(E.f,Cr,N,wd.m, wl.m, q.m,efs.m,alpha.m,beta.m, ret.m, rho, restriction){
        # if catch = TRUE (=> the restriction is catch not landings. )

        Ca.m <- CobbDouglasAge(E = E.f, N = N, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m, q.m = q.m,
                                 efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, rho = rho)

        if(restriction == 'catch') Ca.m <- Ca.m
        else  Ca.m <- ret.m*Ca.m
        
                return(Cr - sum(Ca.m))
    }


    Cinf <- CobbDouglasAge(E = 1e100,Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m,  ret.m = ret.m, rho = rho)
    if((Cr - sum(Cinf))> 0) # Even with infinity effort it is not possible to catch the quota => return 'almost' infinity effort.
        return(effort = 1e100)

    NomEff <- uniroot(fObj,interval=c(0,1e100),Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, rho = rho,  restriction = restriction, ret.m = ret.m)$root

    return(effort =  NomEff)
}