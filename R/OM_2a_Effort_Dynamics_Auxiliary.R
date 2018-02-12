#------------------------------------------------------------------------------------
# Dorleta Garcia
# 2017/02/05
#
# Transform the FLR objects into list of arrays in order to be able to work with non-FLR
# R functions.
#------------------------------------------------------------------------------------


FLObjs2S3_fleetSTD <- function(biols, fleets, BDs, advice, covars, biols.ctrl, fleets.ctrl,  flnm, yr, ss,iters){
  
  # Fleet's info
  fl    <- fleets[[flnm]]
  sts   <- catchNames(fl)
  mtnms <- names(fl@metiers)
  nmt   <- length(mtnms)
  nst <- length(biols)
  
  i <- iters
  
  
  
    # Biomass at age.: B[nst,it]
    #------------------------------
    B    <- matrix(t(sapply(names(biols), function(x){   # biomass in the middle of the season  [nst,it]
                                      
                                if(dim(biols[[x]]@n)[1] > 1)
                                  return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss,,i, drop=T])
                                else{
                                  if(biols.ctrl[[x]][['growth.model']] == 'fixedPopulation'){
                                    return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss,,i, drop=T])
                                  }
                                  else{
                                    return((biols[[x]]@n*biols[[x]]@wt + BDs[[x]]@gB)[,yr,,ss,,i, drop=T])
                                  }
              } })) , nst,1, dimnames = list(names(biols), 1:1))

    
    # Numbers at age.: list(nst)[na_st,1,nu_st,1,1,it]
    # biomass at age in the middle  of the season, list elements: 
    #-----------------------------------------------
    N   <- lapply(names(biols), function(x){  
                          if(dim(biols[[x]]@n)[1] > 1)
                              return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss,,i, drop = FALSE])
                          else{
                              if(biols.ctrl[[x]] == 'fixedPopulation'){
                                return((biols[[x]]@n)[,yr,,ss,,i, drop=F])
                              }
                              else{
                                return((biols[[x]]@n + BDs[[x]]@gB)[,yr,,ss,,i, drop=F])
              } } })
    names(N) <- names(biols)

    # Calculate QS by fleet for the year and season: matrix [nf,nst]
    #------------------------------------------------------------------
    QS.fls   <- sapply(names(biols), function(x){           # 
                    yr.share    <- advice$quota.share[[x]][,yr,,,,i, drop=T]        # [nf]
                    ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss,,i, drop=T]   # [nf]
                    quota.share <- yr.share*ss.share # [nf]
                    quota.share[is.na(quota.share)] <- 0
                    return(quota.share)})
    names(QS.fls) <- names(biols)

    # If TAC >= B*alpha => TAC = B*alpha.
    # Total seasonal quota share: [nst] 
    #------------------------------------------------------------------
    TAC.yr   <- advice$TAC[,yr,,,,i,drop=T]    # nst
    rho      <- fleets.ctrl$catch.threshold[,yr,,ss,,i, drop=T]  # [ns]
    QS.ss    <- colSums(QS.fls)  #

    TAC <- ifelse(B*rho < TAC.yr*QS.ss, B*rho, TAC.yr*QS.ss)

    # Re-scale QS to fleet share within the season instead of season-fleet share within year.
    QS   <- sweep(QS.fls, 2, apply(QS.fls,2, sum), "/")  # [nf,nst]
    QS[is.na(QS)] <- 0
    
    # flinfo: matrix with information on which metier catch which stock.
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')
    
    # previous year effort share because we don't know the effort share this yea: [nmt]
    efs.m <- sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr-1,,ss,,i, drop=T]) 
    
    # variable cost per unit of effort [nmt] 
    vc.m <- sapply(mtnms, function(x) fl@metiers[[x]]@vcost[,yr,,ss,,i, drop=T])  #[nmt]
    
    # fixed cost per vessel
    fc    <- fl@fcost[,yr,,ss,,i, drop=T]*covars$NumbVessels[flnm,yr,,ss,,i, drop=T] # [1]
    
    # Crew-share [it]
    crewS <- fl@crewshare[,yr,,ss,,i, drop=T] # [i]

    effs <- numeric(length(sts)); names(effs) <- sts
    Cr.f <- numeric(nst); names(Cr.f) <- names(biols)

    q.m <- alpha.m <- beta.m  <- pr.m <- ret.m <- wd.m <- wl.m <-vector('list', length(sts))
    names(q.m) <- names(pr.m) <- names(alpha.m) <- names(beta.m) <- names(ret.m) <- names(wl.m) <- names(wd.m) <- sts

    tacos <- logical(length(sts))
    names(tacos) <- sts
    
    K <- c(fl@capacity[,yr,,ss,,i,drop=T])

    for(st in sts){     # q.m, alpha.m.... by metier but stock specific
  
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
    Cr.f[st] <- TAC[st,i]*QS[flnm,st]
    }
    
    return(list(B = B,
                N = N,
                QS.fls = QS.fls,
                TAC = TAC,
                rho = rho,
                QS = QS,
                efs.m = efs.m,
                vc.m = vc.m,
                fc = fc,
                crewS = crewS,
                effs = effs, 
                Cr.f = Cr.f,
                tacos = tacos,
                q.m = q.m,
                alpha.m = alpha.m,
                beta.m = beta.m,
                pr.m = pr.m,
                ret.m = ret.m,
                wd.m = wd.m,
                wl.m = wl.m,
                K = K))
}
