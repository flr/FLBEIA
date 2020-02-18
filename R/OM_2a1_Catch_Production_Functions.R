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
#   - Baranov - Catch from Baranov production function by age.
#   - Baranov.effort - Effort from Baranov production.
# 
# Dorleta Garcia 
# Agurtzane Urtizberea (Baranov, January 2020)
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
CobbDouglasBio.effort   <- function(Cr,N, wl.m, wd.m,q.m,efs.m,alpha.m,beta.m,ret.m, rho = NULL, restriction = 'catch',stknm,...){
  
    # if(is.null(rho)){ 
    #   rho <- rep(1, length(N)) 
    #   names(rho) <- names(N)
    # }

  if(is.null(rho))
    rho <- matrix(1, length(N), 1, dimnames = list(names(N), 1))
  
    Cr      <- Cr[stknm,]
    N       <- N[[stknm]] 
    wl.m    <- wl.m[[stknm]]
    wd.m    <- wd.m[[stknm]]
    ret.m   <- ret.m[[stknm]] 
    q.m     <- q.m[[stknm]]
    alpha.m <- alpha.m[[stknm]]
    beta.m  <- beta.m[[stknm]] 
    rho     <- rho[stknm,]

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


    NomEff <- uniroot(fObj,interval=c(0,upl),Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m,  q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, restriction = restriction, ret.m = ret.m, rho = rho,  tol = 1e-12)$root

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

#------------------------------------------------------ -------------------------
#  CobbDouglasAge.Effort :: Cr[1], B[na,nu], efs.m[mt], q.m,alpha.m,beta.m :: [mt,na,nu] 
#-------------------------------------------------------------------------------

CobbDouglasAge.effort   <- function(Cr,N,wl.m, wd.m, ret.m, q.m,efs.m,alpha.m,beta.m, rho = NULL, restriction = 'catch',stknm,...){
   
  
  # if(is.null(rho)){ 
  #   rho <- rep(1, length(N)) 
  #   names(rho) <- names(N)
  # }
  
   if(is.null(rho))
    rho <- matrix(1, length(N), 1, dimnames = list(names(N), 1))
  
     Cr      <- Cr[stknm,]
     N       <- N[[stknm]] 
     wl.m    <- wl.m[[stknm]]
     wd.m    <- wd.m[[stknm]]
     ret.m   <- ret.m[[stknm]] 
     q.m     <- q.m[[stknm]]
     alpha.m <- alpha.m[[stknm]]
     beta.m  <- beta.m[[stknm]] 
     rho     <- rho[stknm,]


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

    NomEff <- uniroot(fObj,interval=c(0,1e100),Cr=Cr,N=N, wl.m = wl.m, wd.m = wd.m, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, rho = rho,  restriction = restriction, ret.m = ret.m, tol = 1e-12)$root

    return(effort =  NomEff)
}



#-------------------------------------------------------------------------------
#  CobbDouglasComb :: 
# If Age structure calls CobbDouglasAge
# If bio structure calls CobbDouglasBio
#-------------------------------------------------------------------------------

CobbDouglasComb   <- function(E,N, wl.m, wd.m, ret.m,q.m,efs.m,alpha.m,beta.m,rho = 1,...){
  
  if(dim(N)[1] == 1) res <- CobbDouglasBio(E = E,N = N, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m,q.m = q.m, efs.m = efs.m,
                                           alpha.m = alpha.m, beta.m  = beta.m, rho = rho)
  else res <- CobbDouglasAge(E = E,N = N, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m,q.m = q.m, efs.m = efs.m,
                             alpha.m = alpha.m, beta.m  = beta.m, rho = rho)

    
  return(res)

}

#-------------------------------------------------------------------------------
#  CobbDouglasAge.Effort :: Cr[1], B[na,nu], efs.m[mt], q.m,alpha.m,beta.m :: [mt,na,nu] 
#-------------------------------------------------------------------------------

CobbDouglasComb.effort   <- function(Cr,N,wl.m, wd.m, ret.m, q.m,efs.m,alpha.m,beta.m, rho = 1, restriction = 'catch', QS.groups, stknm,...){
  
  
  # Identify the stocks in the group of stknm
  grp <- names(which(QS.groups == QS.groups[stknm]))
  
  # The TAC quota of this stocks
  Crs <- sum(unlist(Cr[grp,]))
  
  fObj <- function(E.f,Crs,N,wd.m, wl.m, q.m,efs.m,alpha.m,beta.m, ret.m, rho, restriction, grp = grp){
    # if catch = TRUE (=> the restriction is catch not landings. )
    
    lCa.m <- sapply(grp, function(x){ 
      
            Ca.m <- CobbDouglasComb(E = E.f, N = N[[x]], wl.m = wl.m[[x]], wd.m = wd.m[[x]], 
                                           ret.m = ret.m[[x]], q.m = q.m[[x]], efs.m = efs.m, 
                                           alpha.m = alpha.m[[x]], beta.m = beta.m[[x]], rho = rho[x,])
            if(restriction == 'catch') Ca.m <- Ca.m
            else  Ca.m <- ret.m[[x]]*Ca.m
            
            res <- sum(Ca.m)
            return(res)
    })
    
    
    return(sum(Crs) - sum(lCa.m))
  }
  
  
  Cinfs <- sapply(grp, function(x) 
                      Cinf.x <- sum(CobbDouglasComb(E = 1e100, Cr = Cr[x,], N = N[[x]], wl.m = wl.m[[x]], wd.m = wd.m[[x]], q.m=q.m[[x]],
                                      efs.m=efs.m,alpha.m=alpha.m[[x]],beta.m=beta.m[[x]],  ret.m = ret.m[[x]], rho = rho[x,])))
    
  if((sum(Crs) - sum(Cinfs))> 0) # Even with infinity effort it is not possible to catch the quota => return 'almost' infinity effort.
      return(effort = 1e100)
  
  NomEff <- uniroot(fObj,interval=c(0,1e100),Cr=Crs,N=N, wl.m = wl.m, wd.m = wd.m, q.m=q.m,efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m, rho = rho,  restriction = restriction, ret.m = ret.m, tol = 1e-12, grp = grp)$root
  
  return(effort =  NomEff)
}


#-------------------------------------------------------------------------------
#  Baranov :: E[it], B[na,nu,it], efs.m[mt,it], q.m,alpha.m,beta.m :: [mt,na,nu,it] 
# Res => C.m[mt,na,nu,it=1]
# OUTPUT: Catch at age in weight.
#-------------------------------------------------------------------------------

Baranov   <- function(E, Cr, efs.m, wl.m, wd.m, ret.m, q.m,
                       alpha.m, rho = 1, tac, Cyr_1 = NULL, Nyr_1 = NULL, Myr_1 = NULL, N, M, Cfyr_1 = NULL, Fknown = FALSE, Ft = NULL, Ff= NULL, ...){
  
  
  na <- dim(M)[1]
  W    <- ret.m*wl.m + (1-ret.m)*wd.m
  
  #We estimate F total assuming being proportional to F yr-1 and gamma an escalar
  if(Fknown == FALSE){
    # We apply baranov at fleet level so we take the mean of the weigths.
    Wall <- array(apply(ret.m*wl.m + (1-ret.m)*wd.m,2:4,mean), dim = dim(M))

    findF <- function(Fa,C, M, N){
      C. <- (Fa/(Fa+M))*N*(1-exp(-M-Fa))
      res <- sum(C.)-sum(C)
      return(res)
    }
  
    Fayr_1 <- Nyr_1
    Fayr_1[] <- NA
  
    for(a in 1:na) Fayr_1[a,,] <- uniroot(findF,interval=c(0,2), C = Cyr_1[a,,], M = Myr_1[a,,], N = Nyr_1[a,,],  tol = 1e-12,extendInt = "yes")$root
  
    # estimate gamma where Fa =Fayr-1 *gamma
    #if(na>1) Na  <-  as.vector(N)*iter(exp((biols[[stknm]]@m/2)[,yr,,ss]),it) # N is in the middle of the season,
  
    find_gamma <- function(gamma, Fayr_1,N,M,W,tac){
      z   <- Fayr_1*gamma + M
      C.  <- (Fayr_1*gamma/z) *(1-exp(-z))*N*W
      res <- sum(C.)-tac
      return(res)}
  
  
      gamma <- uniroot(find_gamma, interval=c(-5,5), Fayr_1 = Fayr_1, N = N, M = M, W = Wall, tac = tac,
                   tol = 1e-12, extendInt = "yes")$root
  
      Fafyr_1   <- Nyr_1
      Fafyr_1[] <- NA
  
      for (a in 1:na) Fafyr_1[a,,] <- uniroot(findF,interval=c(0,2), C = Cfyr_1[a,,], M = Myr_1[a,,], N = Nyr_1[a,,], 
                                             tol = 1e-12, extendInt = "yes")$root
  
      find_delta <- function(delta,gamma, Fayr_1, Fafyr_1,N,M,W,Cr){
        z   <- Fayr_1*gamma+M
        Cf  <- delta*Fafyr_1/z*(1-exp(-z))*N*W
        res <- sum(Cf)-Cr
      }
  
      delta<- uniroot(find_delta,interval=c(-5,5), gamma = gamma, Fayr_1 = Fayr_1, Fafyr_1 = Fafyr_1, N = N, M = M, W = Wall, Cr = Cr,
                  tol = 1e-12,extendInt = "yes")$root
  
      #Now Ff <- Ff*delta is known
      dimq <- dim(q.m)
  
      Ef  <- matrix(E,dim(efs.m)[1],dim(efs.m)[2], byrow = T)*efs.m      # [mt,it]
      Ef  <- array(Ef, dim = c(dim(Ef), dimq[2:3]))
      Ef  <- aperm(Ef, c(1,3:4,2))   # [mt,na,nu,it] 
  
      Faf_new <- array(apply(q.m*Ef,2:4, sum), dim = dim(q.m)[-1])
  
      Fa_new <- gamma*Fayr_1-delta*Fafyr_1+ Faf_new
      z_new <- Fa_new+M
      z_new <- array(z_new,dim=c(dim(Ef)))
      N <- array(N,dim=c(dim(Ef)))

      C.m <- q.m*Ef/z_new*(1-exp(-z_new))*N*W
  
      Clim <- sweep(W*N,4,rho,"*")
  
      C.m <- ifelse(C.m < Clim, C.m, Clim) # The truncation  of CobDoug is applied at metier level.
  
      return(catch =  C.m)}  # [mt,na,nu,it]
  
      else{ # We know the effort of all the fleets so we can calculate the real F. 
        C.m <- Ff
        C.m[] <- NA
        for(mt in 1: dim(Ff)[1]) C.m[mt,,,] <- array(Ff[mt,,,], dim = dim(M))/(M+Ft)*(1-exp(-M-Ft))*N*array(W[mt,,,], dim = dim(M))  # [mt,na,nu,it]
        return(C.m)
      }  
}

#------------------------------------------------------ -------------------------
#  Baranov.Effort :: Cr[1], B[na,nu], efs.m[mt], q.m,alpha.m,beta.m :: [mt,na,nu] 
#-------------------------------------------------------------------------------

Baranov.effort   <- function(Cr,wl.m, wd.m, ret.m, q.m,efs.m,alpha.m,beta.m, rho = NULL,
                              restriction = 'catch',stknm,QS.groups,tac,Cyr_1,Nyr_1,Myr_1,N,M,Cfyr_1,...){
  
  if(is.null(rho)){ 
    rho <- rep(1, length(N)) 
    names(rho) <- names(N)
  }
  
  Cr      <- Cr[stknm,]
  N       <- N[[stknm]] 
  wl.m    <- wl.m[[stknm]]
  wd.m    <- wd.m[[stknm]]
  ret.m   <- ret.m[[stknm]] 
  q.m     <- q.m[[stknm]]
  alpha.m <- alpha.m[[stknm]]
  beta.m  <- beta.m[[stknm]] 
  rho     <- rho[stknm,]
  tac     <- tac[stknm,]
  Cyr_1   <- Cyr_1[[stknm]] 
  Nyr_1   <- Nyr_1[[stknm]]    
  Myr_1   <- Myr_1[[stknm]] 
  M       <- M[[stknm]] 
  Cfyr_1  <- Cfyr_1[[stknm]] 
  
  
  fObj <- function(E.f,Cr,wd.m, wl.m, q.m,efs.m,alpha.m,beta.m, ret.m, rho, restriction,
                   tac, Cyr_1=Cyr_1,Nyr_1=Nyr_1,Myr_1=Myr_1,N=N,M=M,W=W,Cfyr_1=Cfyr_1){
    # if catch = TRUE (=> the restriction is catch not landings. )
    
    Ca.m <- Baranov(E = E.f, Cr=Cr, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m, q.m = q.m,
                     efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, rho = rho,
                     tac = tac, Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1, 
                     N = N, M = M, Cfyr_1 = Cfyr_1)
    
    if(restriction == 'catch') Ca.m <- Ca.m
    else  Ca.m <- ret.m*Ca.m
    
    return(Cr - sum(Ca.m))
  }
  
  
  Cinf <- Baranov(E = 1e100,Cr=Cr,wl.m = wl.m, wd.m = wd.m, q.m=q.m,
                   efs.m=efs.m,alpha.m=alpha.m,beta.m=beta.m,  ret.m = ret.m, rho = rho,
                   tac=tac, Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1, N = N, M = M,
                    Cfyr_1 = Cfyr_1)
  if((Cr - sum(Cinf))> 0) # Even with infinity effort it is not possible to catch the quota => return 'almost' infinity effort.
    return(effort = 1e100)
  
  NomEff <- uniroot(fObj,interval=c(0,1e50),Cr=Cr, wl.m = wl.m, wd.m = wd.m, q.m=q.m,
                    efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, rho = rho,
                    restriction = restriction, ret.m = ret.m,
                    tac=tac, Cyr_1 = Cyr_1, Nyr_1 = Nyr_1, Myr_1 = Myr_1, N = N, M = M,
                    W = W, Cfyr_1 = Cfyr_1, tol = 1e-4)$root
  
  return(effort =  NomEff)
}

