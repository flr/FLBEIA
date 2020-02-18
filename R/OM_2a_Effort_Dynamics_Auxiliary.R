#------------------------------------------------------------------------------------
# Dorleta Garcia
# 2017/02/05
# 2019/02/17
#
# Transform the FLR objects into list of arrays in order to be able to work with non-FLR
# R functions.
# Need to check rho and tac and QS.
# rho: Check what sonia has done the last days.
#------------------------------------------------------------------------------------


FLObjs2S3_fleetSTD <- function(biols, fleets, BDs, advice, covars, biols.ctrl, fleets.ctrl,  flnm, yr, ss,iters){
  
  # Fleet's info
  fl    <- fleets[[flnm]]
  sts   <- catchNames(fl)
  stnms <- names(biols)
  flnms <- names(fleets)
  mtnms <- names(fl@metiers)
  nmt   <- length(mtnms)
  nst <- length(biols)
  nsts <- length(sts)
  nit <- length(iters)
  
  i <- iters
  
  effort.model <- fleets.ctrl[[flnm]][['effort.model']]

    # Biomass at age.: B[nst,it]
    #------------------------------
    B    <- matrix(t(sapply(sts, function(x){   # biomass in the middle of the season  [nst,it]
                                      catch.model <- fleets.ctrl[[flnm]][[x]][['catch.model']]
                                      if(catch.model == 'CobbDouglasAge') return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss,,iters, drop=T])
                                      if(biols.ctrl[[x]][['growth.model']] == 'fixedPopulation') return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss,,iters, drop=T])
                                      if(catch.model == 'Baranov') return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt))[,yr,,ss,,iters, drop=T])
                                      return((biols[[x]]@n*biols[[x]]@wt + BDs[[x]]@gB)[,yr,,ss,,iters, drop=T])
                                      })), nsts,nit, dimnames = list(sts, iters))
  
    
    # Numbers at age.: list(nst)[na_st,1,nu_st,1,1,it]
    # biomass at age in the middle  of the season, list elements: 
    #-----------------------------------------------
    N   <- lapply(setNames(sts, sts), function(x){   # biomass at age in the middle  of the season, list elements: [na,1,nu,1,1,it]
                                catch.model <- fleets.ctrl[[flnm]][[x]][['catch.model']]
                                if(catch.model == 'CobbDouglasAge')        return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss,,iters, drop = FALSE])
                                if(biols.ctrl[[x]] == 'fixedPopulation')   return((biols[[x]]@n)[,yr,,ss,,iters, drop=F])
                                if(catch.model == 'Baranov')               return(biols[[x]]@n[,yr,,ss,,iters, drop = FALSE])
                                # else   
                                return((biols[[x]]@n + BDs[[x]]@gB)[,yr,,ss,,iters, drop=F])})

    names(N) <- sts


    # Calculate QS by fleet for the year and season:
    # If maxprofit dimension: matrix [nf,nst]
    # If SMFB: list[nf,it]
    #------------------------------------------------------------------
    if(effort.model == 'MaxProfit'){
      QS   <- sapply(names(biols), function(x){           # 
                    yr.share    <- advice$quota.share[[x]][,yr,,,,iters, drop=T]        # [nf, it] or [nf]
                    ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss,,iters, drop=T]   # [nf, it] or [nf]
                    quota.share <- yr.share*ss.share # [nf, it] or [nf]
                    quota.share[is.na(quota.share)] <- 0
                    return(quota.share)})
      names(QS) <- names(biols)
    }
    else{
      QS   <- lapply(stnms, function(x){           # list of stocks, each stock [nf,it]
        # Calculate QS by fleet for the year and season
        yr.share    <- advice$quota.share[[x]][,yr,, drop=T]        # [nf,it]
        ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss, drop=T]   # [nf,it]
        quota.share <- matrix(yr.share*ss.share, length(flnms), nit, dimnames = list(flnms, iters))
        quota.share[is.na(quota.share)] <- 0
        return(quota.share)})         
        names(QS) <- stnms
    }
    
    
    # The following objects are used in  Baranov. Baranov needs the biomass at the beginning of the season and we need to 
    # make some asumptions about fleet dynamics.
    # All have dimension: [na,1,nu,1,1,it]
    #-----------------------------------------------------
    Nyr_1   <- lapply(setNames(sts, sts), function(x){   # biomass at age beggining of the season yr-1, list elements: [na,1,nu,1,1,it]
      return(biols[[x]]@n[,yr-1,,ss,iters, drop = FALSE])})
    M   <- lapply(setNames(sts, sts), function(x){   # M  yr, list elements: [na,1,nu,1,1,it
      return(biols[[x]]@m[,yr,,ss,,iters,  drop = FALSE])})
    Myr_1   <- lapply(setNames(sts, sts), function(x){   # M yr-1, list elements: [na,1,nu,1,1,it]
      return(biols[[x]]@m[,yr-1,,ss,,iters,  drop = FALSE])})
    Cyr_1  <- lapply(setNames(sts, sts), function(x) catchStock(fleets,x)[,yr-1,,ss,,iters,  drop = FALSE])
    Cfyr_1 <- lapply(setNames(sts, sts), function(x) catchStock.f(fl,x)[,yr-1,,ss,,iters, drop = FALSE])
    
    
    # If TAC >= B*alpha => TAC = B*alpha.
    # if rho is a numeric => there is only one stock and one iteration wihout names => name it
    # [nst,it]
    #------------------------------------------------------------------------------------------
    TAC.yr  <- matrix(advice$TAC[stnms,yr,drop=T], nst, nit, dimnames = list(stnms, iters))   # [nst,it]
    rho       <- fleets.ctrl$catch.threshold[,yr,,ss,,iters,  drop=T]  # [ns,it]
    if(is.null(dim(rho)) & is.null(names(rho))) names(rho) <- stnms
    if(effort.model == 'MaxProfit')
      QS.ss    <- matrix(colSums(QS), dim(QS)[2], nit, dimnames = list(colnames(QS), iters))  # [nst,it] 
    else 
      QS.ss   <- matrix(t(sapply(stnms, function(x) apply(QS[[x]],2,sum))), nst, nit, dimnames = list(stnms, iters))  # [nst,it]
    
    # TAC [nstnms, it], limited by rho*B
    #---------------------------------------------------------
    if(nit > 1){    
      if(length(stnms) == 1) rho <- matrix(rho, 1,nit, dimnames = list(stnms, 1:iters))
      
      TAC <- ifelse(B[sts,]*rho[sts,] < TAC.yr[sts,]*QS.ss[sts,], B[sts,]*rho[sts,], TAC.yr[sts,]*QS.ss[sts,])
    }
    else TAC <- ifelse(B[sts,]*rho[sts] < TAC.yr[sts,]*QS.ss[sts,], B[sts,]*rho[sts], TAC.yr[sts,]*QS.ss[sts,])
    
    TAC <- matrix(TAC, length(sts),nit,dimnames = list(sts, iters))
    

    # Re-scale QS to fleet share within the season instead of season-fleet share within year.
    #------------------------------------------------------------------------------------------
    if(effort.model == 'MaxProfit'){ 
      QS   <- sweep(QS, 2, apply(QS,2, sum), "/")  # [nf,nst]
      QS[is.na(QS)] <- 0
    }
    else{
      QS   <- lapply(stnms, function(x){          # list of stocks, each stock [nf,it]
        res <- sweep(QS[[x]], 2, apply(QS[[x]],2, sum), "/")
        res[is.na(res)] <- 0 
        return(res)})      
      names(QS) <- stnms
    }
    
    
    # The effort is restricted only by the stocks in 'stocks.restr'    
    # If the restrictors are missing => all the stocks restrict.
    #-----------------------------------------------------------------------------------
    if(is.null(fleets.ctrl[[flnm]][['stocks.restr']]) |  length(fleets.ctrl[[flnm]][['stocks.restr']]) == 0) {
      fleets.ctrl[[flnm]][['stocks.restr']] <- catchNames(fleets[[flnm]])
    }  
    sts <- fleets.ctrl[[flnm]][['stocks.restr']]
    
    
    
    # flinfo: matrix with information on which metier catch which stock.
    #-------------------------------------------------------------------------------------
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')
    
    
    # effort share: In MaxProift previous year effort share because we don't know the effort share this year: [nmt]
    #-----------------------------------------------------------------------------------------
    yr.efs <- ifelse(effort.model == 'MaxProfit', yr-1, yr)
    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr.efs,,ss,,iters, drop=T])), 
                            length(mtnms), nit, dimnames = list(metier = mtnms, iters))
    
    
    # Economic input factors to compute profits
    #-----------------------------------------------------------------------------------------
    vc.m <- fc <- crewS <- NULL
    if(effort.model == 'MaProfit'){
      # variable cost per unit of effort [nmt] 
      vc.m <- sapply(mtnms, function(x) fl@metiers[[x]]@vcost[,yr,,ss,,iters, drop=T])  #[nmt]
    
      # fixed cost per vessel
      fc    <- fl@fcost[,yr,,ss,,i, drop=T]*covars$NumbVessels[flnm,yr,,ss,,iters, drop=T] # [1]
    
      # Crew-share [it]
      crewS <- fl@crewshare[,yr,,ss,,iters, drop=T] # [i]
    }


    effs <- matrix(NA,length(sts), nit, dimnames = list(sts, iters))
    Cr.f <- matrix(NA,length(sts), nit, dimnames = list(sts, iters))

    q.m <- alpha.m <- beta.m  <- pr.m <- ret.m <- wd.m <- wl.m <-vector('list', length(sts))
    names(q.m) <- names(pr.m) <- names(alpha.m) <- names(beta.m) <- names(ret.m) <- names(wl.m) <- names(wd.m) <- sts

    tacos <- logical(length(sts))
    names(tacos) <- sts
    
    K <- c(fl@capacity[,yr,,ss,,iters,drop=T])
    
    # q.m, alpha.m.... by metier but stock specific
    #----------------------------------------------------------------------
    for(st in sts){    
  
      # identify the first metier that catch stock st
      mtst <- flinfo[[st]][2]
  
      age.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[1]]
      age.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[1]]
      age.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[1]]
      age.pr    <- dimnames(fl@metiers[[mtst]]@catches[[st]]@price)[[1]]
  
      unit.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[3]]
      unit.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[3]]
      unit.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[3]]
      unit.pr    <- dimnames(fl@metiers[[mtst]]@catches[[st]]@price)[[3]]
  
      q.m[[st]]     <- array(0, dim = c(length(mtnms), length(age.q),     length(unit.q), 1),      dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1))
      alpha.m[[st]] <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), 1), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1))
      beta.m[[st]]  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
      ret.m[[st]]   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
      wl.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
      wd.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
      pr.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), 1),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta, iter = 1))
  
  
      # if TAC overshoot is not discarded, it is sold and it contributes to the revenue.
      tacos[st] <- ifelse(is.null(fleets.ctrl[[flnm]][[st]][['discard.TAC.OS']]), TRUE,fleets.ctrl[[flnm]][[st]][['discard.TAC.OS']])
  
      for(mt in mtnms){
        if(!(st %in% names(fl@metiers[[mt]]@catches))) next
    
        q.m[[st]][mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, ,i,drop = TRUE]
        alpha.m[[st]][mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, ,i,drop = TRUE]
        beta.m[[st]][mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, ,i,drop = TRUE]
        ret.m[[st]][mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, ,i,drop = TRUE]
        wl.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, ,i,drop = TRUE]
        wd.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, ,i,drop = TRUE]
    
        # The price is taken from the year before, because price for the current year is updated after catch is produced,
        # if the price was dynamically updated inside this function the optimizer could crash.
        pr.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@price[,yr-1,,ss,,i, drop = TRUE]
      }
      
      
      # identify the first metier that catch stock st
      mtst <- flinfo[[st]][2]
      
      if(effort.model == 'MaxProfit') Cr.f[st,] <- TAC[st,]*QS[flnm,st]
      else Cr.f[st,] <- TAC[st,]*QS[[st]][flnm,]
      
      LO <- ifelse(length(fleets.ctrl[[flnm]]$LandObl)==1, fleets.ctrl[[flnm]]$LandObl, fleets.ctrl[[flnm]]$LandObl[yr])
      
      if(LO){
        if(fleets.ctrl[[flnm]]$LandObl_yearTransfer[yr-1] == TRUE){ # If landing Obligation = TRUE discount possible quota transfer from previous year.
          Cr.f[st,] <- Cr.f[st,] - fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[st,yr-1,]
          Cr.f[st,] <- ifelse(Cr.f[st,]<0, 0, Cr.f[st,])  # if lower than 0 , set it to 0.
        }
      }

    }
    
    return(list(B = B,
                N = N,
                QS = QS,
                TAC = TAC,
                rho = rho,
                efs.m = efs.m,
                vc.m = vc.m,
                fc = fc,
                crewS = crewS,
                effs = effs, 
                Cr.f = Cr.f,
                TAC.yr = TAC.yr,
                tacos = tacos,
                q.m = q.m,
                alpha.m = alpha.m,
                beta.m = beta.m,
                pr.m = pr.m,
                ret.m = ret.m,
                wd.m = wd.m,
                wl.m = wl.m,
                K = K,
                Nyr_1 = Nyr_1,
                Myr_1 = Myr_1,
                M = M,
                Cfyr_1 = Cfyr_1,
                Cyr_1  = Cyr_1,
                LO = LO))
}
