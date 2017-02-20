# Changes: line 13 (stksnms). Added: lines 39 and 40 (Itsaso Carmona, 01/04/2015)
# Changes: line 13 removed stknms, fleets is not an argument to the function and 
#           the RCMD Check failed, now stknms an argument to the function. 2017/02/6
QuotaSwap <- function(stknms, E0, Cr.f,Cr.f_exemp, N, B, efs.m, q.m, alpha.m, beta.m, pr.m = NULL,wl.m, wd.m, ret.m, fc = NULL, vc.m = NULL, crewS = NULL,  K,rho, flnm, fleets.ctrl, stks_OF, approach = 'maxprof'){

    if(!(approach %in% c('maxprof', 'fcube'))) stop("Approach in QuotaSwap function must be 'maxprof' or 'fcube'.")
    
    Cr.f.new <- Cr.f_exemp  # Object to store the new quotas.
    quota_swap_p  <- rep(0,length(N))
    quota_swap_st <- rep(NA,length(N))
    names(quota_swap_p) <- names(quota_swap_st) <- names(N)
    
#     stksnms <- names(N)
    stksnms <-stknms
    LO_stk_grp <- fleets.ctrl[[flnm]]$LO_stk_grp
      
  
    
    # Calculate the catch derived from eff_nloptr and see if any stock is restricting the effort, 
    # NO  => leave the loop
    # YES => start the swap of quotas.
    CE <- numeric(length(stksnms))
    Ca_st <- vector('list', length(CE))
    names(Ca_st) <-   names(CE) <- stksnms
    
    nmt <- dim(efs.m)[1]

    for(st in stksnms){ 
      catchFun <- fleets.ctrl[[flnm]][[st]][['catch.model']]
      Nst      <- N[[st]]
      Ca_st[[st]] <-  (eval(call(catchFun, N = Nst, B = B[st], E = E0, efs.m = efs.m, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]], rho = rho[st])))
      
      if(dim(Ca_st[[st]])[2] == 1) Ca_st[[st]] <- array(Ca_st[[st]], dim = c(dim(Ca_st[[st]])[1:2],1,1), dimnames = list(dimnames(Ca_st[[st]])[[1]],1,1,1)) # unit = iter = 1
      
      CE[st] <- sum(Ca_st[[st]])
      # CE is the catch obtained using the initial Effort that includes the minimis and the quota transfer btw years.
      La_st <- Da_st <- Ca_st
      La_st[[st]] <- Ca_st[[st]]*ret.m[[st]]
      Da_st[[st]] <- Ca_st[[st]] - La_st[[st]]
    }    
    
    E <- E0 # just in case the condition in 'while' is false from the start.
    
names(Cr.f)<- stksnms
names(Cr.f.new)<- stksnms

Cr.f.new <- ifelse(Cr.f.new == 0, 1e-6, Cr.f.new)
Cr.f <- ifelse(Cr.f == 0, 1e-6, Cr.f)

    while(any((CE[stksnms]/Cr.f.new[stksnms]) > 0.98) & any((CE[stksnms]/Cr.f[stksnms]) < 0.98)){ # There is al least one stock which is restricting the effort.  OR THE Quota compsumtion of **ALL** the stocks is above 98%
        
        STRs <- stksnms[which((CE[stksnms]/Cr.f.new[stksnms]) > 0.98)] #  RESTRICTORS, we start increasing the quota of the first stock. 
        
        STDs        <- vector('list', length(STRs))
        names(STDs) <- STRs
        
        tau.old <- numeric(length(STRs))
        names(tau.old) <- STRs
        
     #   browser()
        
        lambda.lim <- (1 - CE[stksnms]/Cr.f.new[stksnms]) # The surplus for each of the stocks.
 #  print(lambda.lim)
          
        for(str in STRs){
    
            # Get the group of donors
            if(!is.na(quota_swap_st[str])){ 
                STDs[[str]]    <- quota_swap_st[str]
            }
            else{ # Is the first time that the stocks receives a 'donation'.
                STDs[[str]] <- names(LO_stk_grp[(LO_stk_grp == LO_stk_grp[str]) & !(names(LO_stk_grp) %in% STRs)])
            }
            
            # Maintain in the list only the stocks that:

            #  * Can donate more than a 2% of their **original** quota share (with this restriction we would be excluding the stocks that have received quota from minimise, year transfer or donation)  
            STDs[[str]] <- stksnms[which((CE[STDs[[str]]]/Cr.f[STDs[[str]]]) < 0.98)]
            
            if( (quota_swap_p[str] > 0.089) | all(lambda.lim[STDs[[str]]] < 0.02) | length(STDs[[str]]) == 0 | stks_OF[str] == TRUE) 
              return(list(E = E, efs.m = efs.m, Cr.f.new = Cr.f.new, quota_swap_st = quota_swap_st, quota_swap_p = quota_swap_p, catch = CE, Ca = Ca_st, La = La_st, Da = Da_st))
                 }
 
 print(STDs)

        STDs.grid <- expand.grid(STDs)
            
        std_opt_pars <- matrix(NA,ifelse(approach == 'maxprof', nmt+length(STRs),length(STRs)), dim(STDs.grid)[1]) 
        std_opt_obj <- numeric(dim(STDs.grid)[1])
        
        tau.old <- quota_swap_p[STRs]
      
        for(k in 1:dim(STDs.grid)[1]){   
          
            STDs <- as.matrix(STDs.grid[k,])
            colnames(STDs) <- names(STDs.grid)
              
            print(STRs)
            print(STDs)
       
    #        print(lambda.lim)
            
            tau.init  <- numeric(length(STRs))
            names(tau.init) <- STRs
            tau.up <- tau.init
    
            for(st in unique(STDs[1,])){
              tau.init[names(which(STDs[1,] == st))] <- ifelse((lambda.lim[st]/length((which(STDs[1,] == st)))) >= 0.09, 0.089, lambda.lim[st]/length((which(STDs[1,] == st))))  
              tau.up[names(which(STDs[1,] == st))] <- ifelse(lambda.lim[st] > 0.09, 0.09, lambda.lim[st])  
            }
    
  #  print(tau.up)
          
            if(approach == 'maxprof'){
               E0m <- E0*efs.m
               
        #       browser()
                
                Nnl <- lapply(N, function(x) array(x, dim = c(dim(x)[1],1,dim(x)[2],1,1,dim(x)[3])))
               
                std_nloptr <- nloptr::nloptr(c(tau.init,E0m),
                           eval_f= f_MP_LO_nloptr,
                           lb = c(rep(0,length(STRs)), E0m*0.999),
                           ub = c(tau.up, E0m + (K-E0)),
                          eval_g_ineq = g_ineq_MP_LO_nloptr ,
                           opts = list("algorithm" = "NLOPT_LN_COBYLA", maxeval = 100, xtol_abs = 1e-4, xtol_rel = 1e-4, maxtime = 300),
                           q.m = q.m, alpha.m = alpha.m, beta.m = beta.m, pr.m = pr.m,  Cr.f = Cr.f, fc = fc,
                           ret.m = ret.m, wd.m = wd.m, wl.m = wl.m, vc.m = vc.m, N = Nnl,  B = B,  K=K,  rho = rho,
                           crewS = crewS, STRs = STRs, STDs = STDs,  Cr.f.new = Cr.f.new, tau.old = tau.old, lambda.lim = lambda.lim)
         
            }
            else{ # => approach == fcube

                
              std_nloptr <- nloptr::nloptr(tau.init,
                                   eval_f= f_fcube_LO_nloptr,
                                   lb = rep(0,length(STRs)),
                                   ub = tau.up,
                                   eval_g_ineq = g_ineq_fcube_LO_nloptr ,
                                   opts = list("algorithm" = "NLOPT_LN_COBYLA", maxeval = 100, xtol_abs = 1e-4, xtol_rel = 1e-4, maxtime = 300),
                                   q.m = q.m, alpha.m = alpha.m, beta.m = beta.m,  Cr.f = Cr.f, efs.m = efs.m,
                                   ret.m = ret.m, wd.m = wd.m, wl.m = wl.m,  N = N,  B = B,  K=K,  rho = rho,
                                   STRs = STRs, STDs = STDs,  Cr.f.new = Cr.f.new, tau.old = tau.old, lambda.lim = lambda.lim)
        
            }
    
            std_opt_pars[,k]   <- std_nloptr$solution
            std_opt_obj[k]    <- -std_nloptr$objective
     #       print(round(std_nloptr$solution,6)[1:length(STRs)])      

        }
 
      #  print(round(std_opt_obj,3))
        
 
        pars <- std_opt_pars[,which.max(std_opt_obj)]
        tau.new <- pars[1:length(STRs)]
        names(tau.new) <- STRs
        
        if(approach == 'maxprof'){
            E     <- sum(pars[-(1:length(STRs))])
            efs.m <- matrix(pars[-(1:length(STRs))]/E, ncol = 1)
        }
        else{
            E <- max(std_opt_obj) # efs.m in this case is the input to the QuotaSwap function.
        }
        
        STDs    <- as.matrix(STDs.grid[which.max(std_opt_obj),])
        colnames(STDs) <- names(STDs.grid)
        
        quota_swap_st[STRs] <- STDs[,STRs] 
        quota_swap_p[STRs]  <- tau.new + tau.old 

# print(quota_swap_st)
        # Calculate the catch corresponding with new Effort and update Cr.f.new
        # Divide de catch in 'discards (MLS)' and 'landings (NO MLS)'

 #   browser()
        Ca_st <- vector('list', length(stksnms))
        names(Ca_st) <- stksnms
 
        Da_st <- Ca_st
        La_st <- Ca_st

        for(st in stksnms){ 
    
          catchFun <- fleets.ctrl[[flnm]][[st]][['catch.model']]
          Nst      <- N[[st]]
          Ca_st[[st]] <-  (eval(call(catchFun, N = Nst, B = B[st], E = E, efs.m = efs.m, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]], rho = rho[st])))
          if(dim(Ca_st[[st]])[2] == 1) Ca_st[[st]] <- array(Ca_st[[st]], dim = c(dim(Ca_st[[st]])[1:2],1,1), dimnames = list(dimnames(Ca_st[[st]])[[1]],1,1,1)) # unit = iter = 1

          La_st[[st]] <- Ca_st[[st]]*ret.m[[st]]
          Da_st[[st]] <- Ca_st[[st]] - La_st[[st]]
          
          CE[st] <- sum(Ca_st[[st]])
          
        #  if(st == 'OTH') browser()
          
          
          if(st %in% STRs) Cr.f.new[st]   <- Cr.f.new[st] + tau.new[st]*Cr.f[STDs[,st]] # receptor
          else             Cr.f.new[st]   <- Cr.f.new[st] - sum(tau.new[which(STDs[1,] == st)])*Cr.f[st]
          
        }
      }
  
      
    return(list(E = E, efs.m = efs.m, Cr.f.new = Cr.f.new, quota_swap_st = quota_swap_st, quota_swap_p = quota_swap_p, catch = CE, Ca = Ca_st, La = La_st, Da = Da_st))
}
