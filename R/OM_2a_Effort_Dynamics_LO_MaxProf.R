MaxProfit_Extra_LO <- function(biols, fleets, advice.ctrl, fleets.ctrl, fl, Et.res, efs.res,  efs.min, efs.max, 
                               yr,ss, flnm, it, i, sts, q.m, alpha.m, beta.m, pr.m,  Cr.f, fc, ret.m, wd.m, wl.m, vc.m, N, B, K, rho,
                               effort.restr, crewS, catch.restr, efs.abs, tacos){
  
#if(flnm == 'DTS_SP' & yr == 39) browser()
 
  eff <- numeric(it)
  discount_yrtransfer <- matrix(0,length(sts),it, dimnames = list(sts,1:it))
  ret.m.new <- ret.m # retention may change derived from minimis exemption.
  min_ctrl <- rep(FALSE, length(sts))
  names(min_ctrl) <- sts
  
  stnms <- names(biols)
  nmt <- length(fleets[[flnm]]@metiers)

  # Identify the stocks that are unable to 'donate' due to overfishing.
  stks_OF <- overfishing(biols, fleets, advice.ctrl, yr) # matrix[nst,it]

  if(Et.res[i] > c(fl@capacity[,yr,,ss,,i,drop=T])){ 
    fl@effort[,yr,,ss,,i] <- fl@capacity[,yr,,ss,,i,drop=T] 
    next
  }
  else{ # Minimis, Quota transfer btw years and QuotaSwap.
      minimis <- fleets.ctrl[[flnm]]$LandObl_minimis # logical(ny)
      yrtrans <- fleets.ctrl[[flnm]]$LandObl_yearTransfer # logical(ny)
  
      Cr.f_min_qt <- Cr.f
      Et1.res   <- Et.res[i]
      efs1.res <- efs.res[,i]
      
      if(minimis[yr] == TRUE | yrtrans[yr]  == TRUE){
    
        # Add the minimis and quota.transfer 'extra' quota.
        min_p <- fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr]      # nst
        yrt_p <- fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[,yr] # nst
    
        if(length(min_p)==1) names(min_p) <- sts
        if(length(yrt_p)==1) names(yrt_p) <- sts
    
        Cr.f_min_qt[sts] <- Cr.f[sts]*(1+min_p[sts]+yrt_p[sts]) # The quota restriction is enhanced in the proportion allowed by minimis and year transfer.
    
        
        Cr.f_min_qt[sts]   <- ifelse(Cr.f_min_qt[sts] == 0, 1e-8, Cr.f_min_qt[sts])
        
        efs1.res <- ifelse(efs1.res < efs.min, efs1.res, efs.min*1.01)
        efs1.res <- ifelse(efs1.res > efs.max, efs.max*0.99, efs1.res)
        
        Et1.res <- Et1.res*efs1.res
        
        X <- log((Et1.res*efs1.res)/(K - (Et1.res*efs1.res)))
        
        eff_opt <- optim(X,f_MP_nloptr_penalized, efs.max = efs.max, efs.min = efs.min,q.m = q.m, alpha.m = alpha.m, 
                         beta.m = beta.m, pr.m = pr.m, ret.m = ret.m, wd.m = wd.m,
                         wl.m = wl.m, N = N, B = B, fc = fc, vc.m = vc.m,   Cr.f = Cr.f_min_qt,  crewS = crewS, K = K , 
                         effort.restr = 'min', catch.restr = 'catch', efs.abs = efs.abs, tacos = tacos, rho = rho)
        
        
        res <- K/(1+exp(-eff_opt[[1]]))
        Et1.res   <- sum(res)
        efs1.res <- res/sum(res)
 
        # eff_nloptr <- nloptr::nloptr(Et1.res[i]*efs1.res[,i],
        #                          eval_f= f_MP_nloptr,
        #                          lb = rep(0, nmt),
        #                          ub = rep(K, nmt),
        #                          eval_g_ineq = g_ineq_MP_nloptr,
        #                          opts = opts,
        #                          q.m = q.m, alpha.m = alpha.m, beta.m = beta.m, pr.m = pr.m,  Cr.f = Cr.f_min_qt, fc = fc,
        #                          ret.m = ret.m, wd.m = wd.m, wl.m = wl.m, vc.m = vc.m, N = N,  B = B,  K=K,  rho = rho,
        #                          effort.restr = effort.restr, crewS = crewS, catch.restr = restriction, tacos = tacos)
        # Et1.res[i]   <- sum(eff_nloptr$solution)
        # efs1.res[,i] <- eff_nloptr$solution/sum(eff_nloptr$solution)
    
        # Et1.res The effort resulting from minimis and year quota transfer examptions.
        # We will use this effort later to divide the extra catch, in discards (from minimis), year quota transfer 
        # to discount in the following year and quota swap (in this order)
      }
  
      # Quota Swap
      if(!is.null(dim(rho))) rhoi <- rho[,i]
      else rhoi <- rho
  
      # Convert N to the rigth dimension
      #   browser()
      Nqs <- lapply(N,function(x) return(array(x[,1,,1,1,], dim = dim(x)[c(1,3,6)])))
      MP_LO <- QuotaSwap(stknms = sts, E0 = sum(Et1.res), Cr.f = Cr.f, Cr.f_exemp = Cr.f_min_qt, N = Nqs, B = B, efs.m = matrix(efs1.res, nmt), q.m = q.m, 
                     alpha.m = alpha.m, beta.m = beta.m, pr.m = pr.m, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m, 
                     fc = fc, vc.m = vc.m, crewS = crewS, K = K, rho = rho, stks_OF = stks_OF[,i],
                     flnm = flnm, fleets.ctrl = fleets.ctrl, approach = 'maxprof')
  
      efs.res[,i] <- MP_LO$efs.m
      Et.res[i]   <- sum(MP_LO$E)
      cat('LO SWAP: Effort share: ', efs.res[,i], ', ~~~~~ Effort: ',Et.res[i], '\n')
  
      # Divide the extra catch, in discards (from minimis), year quota transfer 
      # to discount in the following year and quota swap (in this order)
      # discount_yrtransfer must be discounted from the quota next year.
  
      catch_Elo <- MP_LO$catch
      diff      <- catch_Elo[sts]/Cr.f[sts] #[nst]
      discount_yrtransfer[sts,i] <- ifelse(diff < 1 + fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr], 0, 
                                       ifelse((diff - fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr] - 1) < fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[,yr], 
                                              (diff - fleets.ctrl[[flnm]]$LandObl_minimis_p[,yr] - 1),
                                              fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[,yr]))*Cr.f[sts]
  
      # update ret.m to account for the discards due to minimise exemption.
      for(st in sts){
    
        #        if(st == 'OTH')
        #          browser()
        # if discards due to size are higher than discards allowed by minimise, ret.m.i is not changed,
        # otherwise is increase so that the total discards equal to min_p*Cr.f  
        Ca <- MP_LO$Ca[[st]]
        Ds <- sum((1-ret.m[[st]])*Ca*wd.m[[st]])                
        ret.m.new[[st]] <- ret.m[[st]] - ifelse(Ds/Cr.f[st] > fleets.ctrl[[flnm]]$LandObl_minimis_p[st,yr], 0, fleets.ctrl[[flnm]]$LandObl_minimis_p[st,yr] - Ds/Cr.f[st])
        min_ctrl[st] <- ifelse(Ds/Cr.f[st]  > fleets.ctrl[[flnm]]$LandObl_minimis_p[st,yr], FALSE, TRUE)
    }
  
  }
  
  # Update the retention curve according to minimis.

  if(any(min_ctrl)){
    sts_min <- names(which(min_ctrl))
  
    for(mt in names(fl@metiers)){
      if(any(sts_min %in% catchNames(fl@metiers[[mt]]))){
        for(st in sts_min[which(sts_min %in% catchNames(fl@metiers[[mt]]))]){
          fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,,,i] <- ret.m.new[[st]][mt,,,]
          fl@metiers[[mt]]@catches[[st]]@discards.sel[,yr,,,,i] <- 1-ret.m.new[[st]][mt,,,]
        }
      }    
    }
  
    fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[,yr,] <- discount_yrtransfer
  }
  
  return(list(fleets.ctrl = fleets.ctrl, fl = fl, discount_yrtransfer = discount_yrtransfer, Et.res = Et.res))
  
}  
  
  