#-------------------------------------------------------------------------------
# Objective function and nonlinear constraint function to calculate the wquota swpas under Fcube like approach.
#
# Dorleta Garcia - Azti Tecnalia
# Created: 22/10/2014 19:02:56
# Changed: 23/10/2014 10:16:01
#-------------------------------------------------------------------------------



#******************************************************************************************
#-------------------------------------------------------------------------------
#  FUNCTION FOR OPTIMIZATION IN LANDING OBLIGATION SCENARIOS
#
#-------------------------------------------------------------------------------
#******************************************************************************************


#-------------------------------------------------------------------------------
# f_LO_nloptr(tau) :: Objective function
#                       (function to be maximized in a fleet by fleet case)
#    - tau: the percentage of quota swaped to restrictive stock.
#    - qa.mt ~ alpha.a.mt ~ beta.a.mt ~ pra.mt.i :: lists with one element per
#                   stock, each element with the form: array(na_st,nmt,it)
#    - Ba ~ list with one element per stock, each element with the form: array(na_st,it)
#
# q.m,...,vc.m must correspond with one iteration of the variable at metier level
# and must have dimension  [nmt,na,nu,1]
# N: alist with one element per stock and each element with dimension [nmt,na,nu]
# Cr.f:[nst,1] quota share
# rho: [nst]
# 
#  Effort restriction is always 'min', i.e all the Quotas must be fulfiled.
# We assuem that the productionfunction is CobbDouglasAge if age structured stock 
# and CobbDouglasBio if bio structured stock
#-------------------------------------------------------------------------------

f_fcube_LO_nloptr <- function(tau, efs.m, q.m, alpha.m, beta.m,  ret.m, wd.m,
                        wl.m, N, B,    Cr.f,  K, rho,STRs,STDs, Cr.f.new, tau.old, lambda.lim){

  names(tau) <- STRs
  
  nmt <- dim(efs.m)[1]
  
  res <- 0
  
  Cst <- Lst <-  numeric(length(q.m))
  names(Cst) <- names(Lst) <-names(q.m)
  
  
  # NEW quotas for str and std
  for(str in STRs){
    Cr.f.new[str] <- Cr.f[str] + tau[str]*Cr.f[STDs[,str]]
    Cr.f.new[STDs[,str]] <- Cr.f.new[STDs[,str]] - (tau[str] - tau.old[str])*Cr.f[STDs[,str]]
 #   cat(str, ' - ', tau,'\n')
  #  print(tau)
  }  
 
  # Calculate the effort corresponding to each stock.
  effs <- numeric(length(N))
  names(effs) <- names(N)
  
  for(st in names(N)){  
     
      Nst <- N[[st]]
      effort.fun <- ifelse(dim(Nst)[1] == 1, 'CobbDouglasBio.effort', 'CobbDouglasAge.effort')
      effs[st] <- eval(call(effort.fun, N = Nst, wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]], Cr = Cr.f.new[st],
                            q.m = q.m[[st]], efs.m = efs.m, alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], rho = rho[st],
                           restriction = 'catch'))

  }    

  E <- min(effs)
  return(-E) # we want to maximize the minimum E.
}

#-------------------------------------------------------------------------------
# Constraints in the swap of quotas in MP.
#   o Total effort in each metier is bigger than the a certain level. (previous effort in the algorithm)
#   o The catch of the stocks != 'str' and 'std' is lower than the 'new catch quota (catch in Cr.f.new).
#   o tau_old < tau_new < 0.9 and:
#       - Cstr < Qstr + (tau_old + tau_new)*Qstd
#       - Cstd < QNEWstd - tau_new*Qstd  
#  Effort restriction is always 'min', i.e all the Quotas must be fulfiled.
#-------------------------------------------------------------------------------
g_ineq_fcube_LO_nloptr <- function(tau, efs.m, q.m, alpha.m, beta.m, ret.m, wd.m, wl.m, N, B, 
                                 Cr.f,  K, rho, STRs,STDs, Cr.f.new, tau.old, lambda.lim){

 # browser()
  
  nmt <- dim(efs.m)[1]
  
  

  E <- - f_fcube_LO_nloptr(tau, efs.m, q.m, alpha.m, beta.m,  ret.m, wd.m,
                                       wl.m, N, B,    Cr.f,  K, rho,STRs,STDs, Cr.f.new, tau.old, lambda.lim)

  stnms <- names(N)
  
  res <- 0
  
  Cst <- Lst <-  NULL
  #   cat( '**************************************************************************\n')
  # constraint on catches, comply with all the TACS ('min') or only with one.
 
  resTAC <- rep(0, length(q.m))  # One resctriction per stock.
  names(resTAC) <- names(q.m)
  
  # NEW quotas for str and std
  for(str in STRs){
      Cr.f.new[str] <- Cr.f[str] + tau[str]*Cr.f[STDs[,str]]
      Cr.f.new[STDs[,str]] <- Cr.f.new[STDs[,str]] - (tau[str] - tau.old[str])*Cr.f[STDs[,str]]
  }  
  
  for(st in names(q.m)){
    Nst <- N[[st]]
    if(dim(Nst)[1] > 1){
        Cam <- CobbDouglasAge(E,Nst, wl.m[[st]], wd.m[[st]],
                              ret.m[[st]],q.m[[st]],matrix(efs.m,ncol = 1),alpha.m[[st]],beta.m[[st]],rho[st])

        resTAC[st] <- sum(Cam)  - Cr.f.new[st]
        #   cat(st, ' - ', sum(Cam), '\n')
      }
    else{
        Cm <- CobbDouglasBio(E,Nst, wl.m[[st]], wd.m[[st]],
                             q.m[[st]],matrix(efs.m,ncol = 1),alpha.m[[st]],beta.m[[st]], ret.m[[st]],rho[st])
        resTAC[st] <- sum(Cm)  - Cr.f.new[st]
        #       cat(st, ' - ', sum(Cm), '\n')
    }

  }
#  browser()
  lambda <- -lambda.lim

  for(st in unique(STDs[1,])){
    strs_std       <- names(STDs[1,which(STDs[1,] == st)])
    lambda[st] <- lambda[st] + sum(tau[strs_std] - tau.old[strs_std])
  }

  resK <- sum(E)-K
# print(c(resTAC, resK, lambda))
# return(c(lambda))
  return(c(resTAC, resK, lambda))
}


