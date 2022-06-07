#-------------------------------------------------------------------------------
#       Update landings, discards, landings.n, discards.n slots.
#
# - aggregated.CobbDoug: Catch at age calculation when Cobb douglas is applied at
#       total biomass level.
# - ageBased.CobbDoug: Catch at age calculation when Cobb douglas is applied at
#       age level.
#  - updateCatch: Update the slots related to catch using the appropiate function.
# Dorleta GarcYYYa
# Created: 28/10/2010 12:33:04
# Changed:03/06/2011 07:53:53
# Changed: 17/02/2020 Baranov added
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# updateCatch(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
updateCatch <- function(fleets, biols, BDs, advice, biols.ctrl, fleets.ctrl, advice.ctrl, year = 1, season = 1){

    fleet.names <- names(fleets)
        
    for(flnm in fleet.names){
        # Which stocks are caught by fleet flnm.
        flsts <- catchNames(fleets[[flnm]])
        for(st in flsts){
            
            na <- dim(biols[[st]]@n)[1]
            if(na == 1) fleets <- BioPop.CAA(fleets = fleets, biols = biols, BDs = BDs, biols.ctrl = biols.ctrl, fleets.ctrl = fleets.ctrl, advice = advice, advice.ctrl = advice.ctrl, year = year, season = season, flnm = flnm, stknm = st)
            else        fleets <- AgePop.CAA(fleets = fleets, biols = biols, BDs = BDs, biols.ctrl = biols.ctrl, fleets.ctrl = fleets.ctrl, advice = advice, advice.ctrl = advice.ctrl, year = year, season = season, flnm = flnm, stknm = st)   
        }
    }
    
    
     fleets <- FLFleetsExt(fleets)
    
    # Correct the catch in case:
    # Age structured: Ca > (Ba*catch.threshold)
    # Biomass: C > (B*catch.threshold)
     
    fleets <- CorrectCatch(fleets = fleets, biols = biols, BDs = BDs, biols.ctrl = biols.ctrl, fleets.ctrl = fleets.ctrl, year = year, season = season)
    
    
    return(fleets)
}




#-------------------------------------------------------------------------------
# aggregated.CobbDoug(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
BioPop.CAA  <- function(fleets, biols, BDs, biols.ctrl, fleets.ctrl, advice, advice.ctrl, year = 1, season = 1, flnm = 1, stknm = 1, ...){

    rho <- fleets.ctrl[['catch.threshold']][stknm,year,,season,drop=T] # [it]

    nf    <- length(fleets)
    stnms <- names(biols)
    nst   <- length(stnms)
    it    <- dim(biols[[stknm]]@n)[6]
    na    <- dim(biols[[stknm]]@n)[1]
    nu <- dim(biols[[stknm]]@n)[3]
    ns    <- dim(biols[[1]]@n)[4]


    if(na > 1) stop('CobbDouglasBio can only be applied at biomass level')

    yr <- year
    ss <- season
    f  <- flnm
    st <- stknm
    
    fleets <- unclass(fleets)

    fl    <- fleets[[f]]
    sts   <- catchNames(fl)
    mtnms <- names(fl@metiers)
    
    adv.ss <- ifelse(is.null(advice.ctrl[[stknm]][["adv.season"]]), ns, advice.ctrl[[stknm]][["adv.season"]])
    
    if(!(st %in% sts)) return(fleets)
    
    # catch restriction, if empty => landings.
    if (is.null(fleets.ctrl[[flnm]]$restriction)) {
      catch.restr <- 'landings'
    } else 
      catch.restr <- ifelse(length(fleets.ctrl[[flnm]]$restriction)==1, fleets.ctrl[[flnm]]$restriction, 
                            fleets.ctrl[[flnm]]$restriction[yr])
    # catch.restr <- ifelse(is.null(fleets.ctrl[[flnm]]$restriction), 'landings',ifelse(length(fleets.ctrl[[flnm]]$restriction)==1, fleets.ctrl[[flnm]]$restriction,fleets.ctrl[[flnm]]$restriction[yr]))
                                             
    tac <- rep(Inf,it)
    
    # if TAC overshoot is discarded, calculate seasonal TAC to calculate the discards.
    TACOS <- fleets.ctrl[[flnm]][[stknm]][['discard.TAC.OS']]
    TACOS <- ifelse(is.null(TACOS), TRUE, TACOS) 
    if(TACOS){
        yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
        ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
        QS          <- yr.share*ss.share                                          # [it]
        QS[is.na(QS)] <- 0
        if (adv.ss < ns & ss <= adv.ss) {
          tac <- advice$TAC[st,yr-1,drop=T]*QS # it
        } else {
          tac <- advice$TAC[st,yr,drop=T]*QS # it
        }
    }
    
    # biomass in the middle of the season  B[it] if age struc.
    if(dim(biols[[st]]@n)[1] > 1){
      B <- unitSums(quantSums(biols[[st]]@n*biols[[st]]@wt*exp(-biols[[st]]@m/2)))[,yr,,ss, drop=T]
      N <- (biols[[stknm]]@n[,yr,,ss]*exp(-biols[[stknm]]@m[,yr,,ss]/2))  # Ba[na,1,1,1,1,it]
    }
    else{ # fixed or  biomass dynamic poulation 
      if(biols.ctrl[[stknm]] == 'fixedPopulation'){
        B <- unitSums(quantSums(biols[[st]]@n*biols[[st]]@wt))[,yr,,ss, drop=T]
        N <- unitSums(quantSums(biols[[st]]@n))[,yr,,ss]
      }
      else{
        B <- BDs[[stknm]]@biomass[,yr,,ss] + BDs[[stknm]]@gB[,yr,,ss]
        N <- B/biols[[stknm]]@wt[,yr,,ss]
      }
      }
    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    eff   <- matrix(fl@effort[,yr,,ss],length(mtnms), it, dimnames = list(mtnms, 1:it), byrow = T)
                     

    # flinfo: matrix with information on which metier catch which stock.
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')

    mtst <- flinfo[[st]][2]

    age.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[1]]
    age.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[1]]
    age.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[1]]

    unit.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[3]]
    unit.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[3]]
    unit.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[3]]

    q.m   <- array(0, dim = c(length(mtnms), length(age.q), length(unit.q),it),     dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1:it))
    alpha.m <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), it), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1:it))
    beta.m  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
    ret.m  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
    wl.m   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
    wd.m   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))


    for(mt in mtnms){

        if(!(st %in% names(fl@metiers[[mt]]@catches))) next

        q.m[mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE]
        alpha.m[mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE]
        beta.m[mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE]
        ret.m[mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, drop = TRUE]
        wl.m[mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, drop = TRUE]
        wd.m[mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, drop = TRUE]
    }
    
    
    Nst  <- array(N[drop=T],dim = dim(N)[c(1,3,6)])
    
    catch.model <- fleets.ctrl[[flnm]][[st]][['catch.model']]  
    
#    cat("Bio - ", flnm, "-", stknm, ": ", catch.model,  "\n")
    
    Cm <- eval(call(catch.model,  E= eff[1,], N = N, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m, q.m = q.m,
             efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, rho = rho))
    
    Ctotal <-  ifelse(rep(catch.restr == 'landings',it), apply(Cm*matrix(ret.m, dim(ret.m)[1], dim(ret.m)[4]),2,sum), apply(Cm,2,sum))

    tac.disc <- ifelse(Ctotal < tac, rep(1,it), tac/Ctotal)
    tac.disc <- ifelse(Ctotal == 0, 0, tac.disc) # In case Ctotal= 0 to avoid NaN in tac.disc

    for(mt in 1:length(mtnms)){

        if(!(st %in% names(fl@metiers[[mt]]@catches))) next

        cobj <- fl@metiers[[mt]]@catches[[st]]

        dsa <- cobj@discards.sel[,yr,,ss]  # [na,1,nu,1,1,it]
        lsa <- cobj@landings.sel[,yr,,ss]  # [na,1,nu,1,1,it]
        sa  <- (dsa + lsa)

        # Recalculate dsa and lsa according to 'tac.disc'     # [na,nu,it]
        lsa <- sweep(lsa,6,tac.disc, "*")
        dsa <- sa - lsa

        cobj@discards[,yr,,ss]   <- Cm[mt,]*dsa # /(sa*tac.disc)
        cobj@landings[,yr,,ss]   <- Cm[mt,]*lsa # *tac.disc/sa
        cobj@discards.n[,yr,,ss] <- cobj@discards[,yr,,ss]/cobj@discards.wt[,yr,,ss]
        cobj@landings.n[,yr,,ss] <- cobj@landings[,yr,,ss]/cobj@landings.wt[,yr,,ss]

        # When land.wt = 0 <-  land.n = NA => change to 0. (idem for disc.wt)
        cobj@landings.n[,yr,,ss][is.na(cobj@landings.n[,yr,,ss])] <- 0
        cobj@discards.n[,yr,,ss][is.na(cobj@discards.n[,yr,,ss])] <- 0

        fl <- fill_flcatches(fl=fl,cobj=cobj,st=st,mt_idx=mt)


        }
          
    fleets[[f]] <- fl
    
#    fleets <- FLFleetsExt(fleets)
      
    return(fleets)
}

            
#-------------------------------------------------------------------------------
# ageBased.CobbDoug(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
AgePop.CAA <- function(fleets, biols, BDs, biols.ctrl, fleets.ctrl, advice, advice.ctrl, year = 1, season = 1, flnm = 1, stknm = 1,...){

    rho <- fleets.ctrl[['catch.threshold']][stknm,year,, season,drop=T] # [it]

    nf    <- length(fleets)
    stnms <- names(biols)
    nst   <- length(stnms)
    it    <- dim(biols[[1]]@n)[6]
    
 #   if(year == 35 & stknm == 'HKE') browser()

    fleets <- unclass(fleets)
    
    yr <- year
    ss <- season
    st <- stknm
    na <- dim(biols[[stknm]]@n)[1]
    nu <- dim(biols[[stknm]]@n)[3]
    ns    <- dim(biols[[1]]@n)[4]
    
    fl    <- fleets[[flnm]]
    sts   <- catchNames(fl)
    mtnms <- names(fl@metiers)
    
    adv.ss <- ifelse(is.null(advice.ctrl[[stknm]][["adv.season"]]), ns, advice.ctrl[[stknm]][["adv.season"]])
    
    catch.model <- fleets.ctrl[[flnm]][[stknm]][['catch.model']]

    if(!(st %in% sts)) return(fleets)
    
    tac <- rep(Inf,it)
 
    # catch restriction, if empty => landings.
    if (is.null(fleets.ctrl[[flnm]]$restriction)) {
      catch.restr <- 'landings'
    } else 
      catch.restr <- ifelse(length(fleets.ctrl[[flnm]]$restriction)==1, fleets.ctrl[[flnm]]$restriction, 
                            fleets.ctrl[[flnm]]$restriction[yr])
    # catch.restr <- ifelse(is.null(fleets.ctrl[[flnm]]$restriction), 'landings',ifelse(length(fleets.ctrl[[flnm]]$restriction)==1, fleets.ctrl[[flnm]]$restriction,fleets.ctrl[[flnm]]$restriction[yr]))
    
 
    #  quota share % to be updated due to year transfer.
    yrtr_p <- fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[stknm,yr]
    yrtr_p <- ifelse(is.null(yrtr_p), 0,yrtr_p)
    # if year transfer was used in previous year discount it, absolute catch
    yrtr_disc <- fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[stknm,yr-1,] # [it]
    
    if(is.null(yrtr_disc)) yrtr_disc <- 0
 
    fleets.ctrl[[flnm]] 
    # if TAC overshoot is discarded, calculate seasonal TAC to calculate the discards.
    TACOS <- fleets.ctrl[[flnm]][[stknm]][['discard.TAC.OS']]    # Is the TAC overshot discarded?
    TACOS <- ifelse(is.null(TACOS), TRUE, TACOS) 
    if(TACOS){
        yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
        ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
        QS          <- yr.share*ss.share                                          # [it]
        QS[is.na(QS)] <- 0              
        if (adv.ss < ns & ss <= adv.ss) {
          tac <- (advice$TAC[st,yr-1,drop=T]*QS*(1+yrtr_p)) - yrtr_disc
        } else {
          tac <- (advice$TAC[st,yr,drop=T]*QS*(1+yrtr_p)) - yrtr_disc   # [it], add yeartransfer in case it is in place, first we increment in % the quota and then we discount the cuota used in previous year. 
                                                                        # the minimise is not added because it is discarded.
        }
    }
    
    if(dim(biols[[st]]@n)[1] == 1) stop(st, ' stock has no ages, Cobb Douglas cannot be applied at age level then! correct the "catch.model" argument in "fleets.ctrl" argument!\n')
    
    if(catch.model == 'Baranov')                  N     <- (biols[[st]]@n[,yr,,ss])  # Ba[na,it], biomass at age in the middle  of the season,
    if(substr(catch.model,1,11) == 'CobbDouglas') N     <- (biols[[st]]@n*exp(-biols[[st]]@m/2))[,yr,,ss]  # Ba[na,it], biomass at age in the middle  of the season,
    
    M     <- biols[[st]]@m[,yr,,ss]  

    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])),
                length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    eff   <- matrix(fl@effort[,yr,,ss],length(mtnms), it, dimnames = list(mtnms, 1:it), byrow = T)
    
    # If 'Baranov' we calculate the total fishing mortality of all the fishery for year 'y'
    if(catch.model == 'Baranov'){
      EFF   <- lapply(setNames(names(fleets), names(fleets)), function(x) c(fleets[[x]]@effort[, yr,,ss]))
      EFS.M <- lapply(setNames(names(fleets), names(fleets)), function(x){ 
        matrix(sapply(names(fleets[[x]]@metiers), function(y)  
                        fleets[[x]]@metiers[[y]]@effshare[, yr,,ss, drop=T] ), length(fleets[[x]]@metiers), it, 
                           dimnames = list(names(fleets[[x]]@metiers), 1:it))})
      
      
       Q.M   <- lapply(setNames(names(fleets), names(fleets)), function(x){ 
                lapply(setNames(names(fleets[[x]]@metiers), names(fleets[[x]]@metiers)),  
                          function(y){
                            if(!(stknm %in% names(fleets[[x]]@metiers[[y]]@catches))) return(array(0, dim = c(na,nu, it)))
                            else return(array(fleets[[x]]@metiers[[y]]@catches[[stknm]]@catch.q[,yr,,ss,drop=TRUE],dim = c(na,nu, it)))
                          })})
       # Calculate total Ft and fleets partial F, Ff.
      FbyFlMt  <- lapply(setNames(names(fleets), names(fleets)), function(x){
                return(lapply(setNames(names(fleets[[x]]@metiers),  names(fleets[[x]]@metiers)), function(y){
                                    Ft.M <- sweep(Q.M[[x]][[y]], 2:3,EFS.M[[x]][y,]*EFF[[x]], "*")
                                    return(Ft.M)}))
        })
      FbyFl  <- lapply(setNames(names(fleets), names(fleets)), function(x){
        return(Reduce( '+', FbyFlMt[[x]]))})
      
      Ft <- array(Reduce('+', FbyFl), dim = c(na, nu, it)) #[na,nu.it]
      Ff <- array(dim = c(length(mtnms), na, nu, it), dimnames = list(mtnms, NULL, 1:nu, 1:it))
      for(m in mtnms) Ff[m,,,] <- FbyFlMt[[flnm]][[m]]  # [mt,na.nu,it]
    }
    else{
      Ft <- Ff <- NULL
    }
    
    # flinfo: matrix with information on which metier catch which stock.
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')

    mtst <- flinfo[[st]][2]

    age.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[1]]
    age.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[1]]
    age.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[1]]

    unit.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[3]]
    unit.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[3]]
    unit.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[3]]

    q.m   <- array(0, dim = c(length(mtnms), length(age.q), length(unit.q),it),     dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1:it))
    alpha.m <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), it), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1:it))
    beta.m  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
    ret.m  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
    wl.m   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
    wd.m   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))


    for(mt in mtnms){

        if(!(st %in% names(fl@metiers[[mt]]@catches))) next

        q.m[mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE]
        alpha.m[mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE]
        beta.m[mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE]
        ret.m[mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, drop = TRUE]
        wl.m[mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, drop = TRUE]
        wd.m[mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, drop = TRUE]
    }

    Nst     <- array(N[drop=T], dim = dim(N)[c(1,3,6)]) # [na.nu.it]
    M       <- array(M[drop=T], dim = dim(M)[c(1,3,6)]) # [na.nu.it]

  #  cat("Age - ", flnm, "-", stknm, ": ", catch.model,  "\n")
    Cam <- eval(call(catch.model, E = eff[1,], N = Nst, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m, q.m = q.m,
                            efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, rho = rho,
                             M = M, Fknown = TRUE, Ft = Ft, Ff = Ff))

    
 # if catch restriction is landings, Lrat is calculated over landigns, else it is calculated over total catch including undersize individuals.
    Ctotal <- ifelse(rep(catch.restr == 'landings', it), apply(Cam*ret.m,4,sum), apply(Cam,4,sum)) 

    tac.disc <- ifelse(Ctotal < tac, rep(1,it), tac/Ctotal)
    tac.disc <- ifelse(Ctotal == 0, 0, tac.disc) # In case Ctotal= 0 to avoid NaN in tac.disc


    Cam <- array(Cam, dim = c(length(mtnms),dim(biols[[st]]@n)[1], 1, dim(biols[[st]]@n)[3],1,it))

    for(mt in 1:length(mtnms)){

        Ca <- array(Cam[mt,,,,,], dim = c(dim(biols[[st]]@n)[1], 1, dim(biols[[st]]@n)[3],1,1,it))

        if(!(st %in% names(fl@metiers[[mt]]@catches))) next
        cobj <- fl[[mt]][[st]]


        na <- dim(q.m)[2]
        nu <- ifelse(is.na(dim(q.m)[3]), 1, dim(q.m)[3])

        efm <- array(eff[1,]*efs.m[mt,], dim = c(it,na,1,nu,1,1))
        efm <- aperm(efm, c(2:6,1))

        dsa <- cobj@discards.sel[,yr,,ss]  # [na,1,nu,1,1,it]
        lsa <- cobj@landings.sel[,yr,,ss]  # [na,1,nu,1,1,it]
        sa  <- (dsa + lsa)  
            
        # Recalculate dsa and lsa according to 'tac.disc'     # [na,nu,it]
        lsa <- sweep(lsa,6,tac.disc, "*")
        dsa <- sa - lsa             # [na,nu,it]

        cobj@discards.n[,yr,,ss] <- Ca*dsa/sa/cobj@discards.wt[,yr,,ss]
        cobj@landings.n[,yr,,ss] <- Ca*lsa/sa/cobj@landings.wt[,yr,,ss]

        # When sa = 0 <-  land.n & dis.n = NA => change to 0.
        cobj@landings.n[,yr,,ss][is.na(cobj@landings.n[,yr,,ss])] <- 0
        cobj@discards.n[,yr,,ss][is.na(cobj@discards.n[,yr,,ss])] <- 0

        cobj@discards[,yr,,ss] <- apply(Ca*dsa/sa,c(2,4,6),sum,na.rm=T)
        cobj@landings[,yr,,ss] <- apply(Ca*lsa/sa,c(2,4,6),sum,na.rm=T)

        fl <- fill_flcatches(fl=fl,cobj=cobj,st=st,mt_idx=mt)      
    }
    
    fleets[[flnm]] <- fl
    
#    fleets <- FLFleetsExt(fleets)
    
    return(fleets)
}


#-------------------------------------------------------------------------------
# seasonshare.CAA(fleets, biols, fleets.ctrl, advice, year = 1, season = 1, flnm = 1, stknm = 1, ...)
#-------------------------------------------------------------------------------
# Estimates catch at age, when catches are derived from fixed season share allocation by metier

seasonshare.CAA  <- function(fleets, biols, fleets.ctrl, advice, advice.ctrl, year = 1, season = 1, flnm = 1, stknm = 1, ...){
  
  # No overshoot allowed
  
  nf    <- length(fleets)
  flnms <- names(fleets)
  stnms <- names(biols)
  nst   <- length(stnms)
  ns    <- dim(biols[[1]]@n)[4] 
  it    <- dim(biols[[1]]@n)[6]
  
  yr <- year
  ss <- season
  f  <- flnm
  st <- stknm
  
  ass.ss <- assess.ctrl[[st]][['ass.season']] #! Need to be corrected as assess.ctrl is not actually an input in the used call
  if (is.null(ass.ss)) { ass.ss <- ns } else if (is.na(ass.ss)) { ass.ss <- ns }
  
  fleets <- unclass(fleets)
  
  fl    <- fleets[[f]]
  sts   <- catchNames(fl)
  mtnms <- names(fl@metiers)
  
  if(!(st %in% sts)) return(fleets)
  
  yy <- ifelse( ss > ass.ss | ass.ss == ns, yr, yr-1)  # for cases when TAC it's not set for natural year
  TAC <- advice$TAC[st,yy]
  
  
  # Find dependencies:
  fl.rel <- fleets.ctrl[[flnm]][[stknm]]$catch.dependence
  if (!is.null(fl.rel)) {
    if (!(fl.rel %in% flnms)) 
      stop("catch.dependence value not valid for '",flnm,"' fleet and '",stknm,"' stock, fleet '",fl.rel,"' not found")
    fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss,] <- fleets.ctrl$seasonal.share[[stknm]][fl.rel,yr,,ss,]
  }
  
  yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
  ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
  QS          <- yr.share*ss.share                                          # [it]
  QS[is.na(QS)] <- 0
  
  Ctotal <- (TAC*QS)[drop=T]
  
  Ba <- biols[[stknm]]@n[,yr,,ss]*biols[[stknm]]@wt[,yr,,ss]*exp(-biols[[stknm]]@m[,yr,,ss]/2)  # Ba[na,1,1,1,1,it]
  B  <- apply(Ba, c(2:6), sum)[drop=T]                                                          # B [it]
  
  CT  <- fleets.ctrl$catch.threshold[st,yr,,ss, drop=T]  # [ns,it]
  
  Ctotal <- ifelse(B*CT < Ctotal, B*CT, Ctotal)
  
  # Check that each metier target only one stock and vice verse
  fl. <- FLFleetsExt(fl); names(fl.) <- flnm
  flinfo     <- stock.fleetInfo(fl.)
  if ( sum(colSums(stock.fleetInfo(fl.))>1)!=0 )
    stop( paste("There is a metier targeting more than one stock, therefore not possible 
                  to use 'seasonshare' catch.model for '",flnm,"' fleet",sep=""))
  if ( sum(rowSums(stock.fleetInfo(fl.))>1)!=0 )
    stop( paste("There is a metier targeting more than one stock, therefore not possible 
                  to use 'seasonshare' catch.model for '",flnm,"' fleet",sep=""))
  flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')
  mt <- flinfo[[st]][2] # metier that captures stock st
  
  cobj <- fl@metiers[[mt]]@catches[[st]]
  
  dsa <- cobj@discards.sel[,yr,,ss]  # [na,1,nu,1,1,it]
  lsa <- cobj@landings.sel[,yr,,ss]  # [na,1,nu,1,1,it]
  sa  <- (dsa + lsa)  
  
  if(dim(biols[[st]]@n)[1] == 1){
    cobj@discards[,yr,,ss]   <- Ctotal*dsa # /(sa*tac.disc)
    cobj@landings[,yr,,ss]   <- Ctotal*lsa # *tac.disc/sa
    cobj@discards.n[,yr,,ss] <- cobj@discards[,yr,,ss]/cobj@discards.wt[,yr,,ss]
    cobj@landings.n[,yr,,ss] <- cobj@landings[,yr,,ss]/cobj@landings.wt[,yr,,ss]
    
    fl@metiers[[mt]]@catches[[st]] <- cobj
  }
  else{ # age structured stock
    #     browser() 
    Bs <- apply(sa*Ba, 6,sum)           # it
    Ca <- sweep(sa*Ba, c(2,4:6), Ctotal/Bs, "*") # [na,nu,it]
    Ca <- ifelse( is.na(Ca), 0, Ca)
    
    cobj@discards.n[,yr,,ss] <- Ca*dsa/sa/cobj@discards.wt[,yr,,ss]
    cobj@landings.n[,yr,,ss] <- Ca*lsa/sa/cobj@landings.wt[,yr,,ss]
    
    # When sa = 0 <-  land.n & dis.n = NA => change to 0.
    cobj@landings.n[,yr,,ss][is.na(cobj@landings.n[,yr,,ss])] <- 0
    cobj@discards.n[,yr,,ss][is.na(cobj@discards.n[,yr,,ss])] <- 0
    
    cobj@landings[,yr,,ss] <- apply(cobj@landings.n[,yr,,ss]*cobj@landings.wt[,yr,,ss],c(2,4,6),sum)
    cobj@discards[,yr,,ss] <- apply(cobj@discards.n[,yr,,ss]*cobj@discards.wt[,yr,,ss],c(2,4,6),sum)
    
    fl@metiers[[mt]]@catches[[st]] <- cobj
  }
  
  fleets[[f]] <- fl
  
  fleets <- FLFleetsExt(fleets)
  
  return(fleets)
}

            
#-------------------------------------------------------------------------------
# CorrectCatch(fleets, biols, year = 1, season = 1)
# Given that in some production functions it can happen that
# Ca > Ba (age struc. pop) or C > B (bio struc. pop), if this happens the
# catch is corrected setting Ca = Ba and C = B. For this end the catch is reduced
# in the same degree in all the fleets. This could imply a 'revision' in the
# production function parameters for the year, season and iteration in question
# but it is not done internally because it does not affect other steps and it
# would imply a loss in the generality of the 'CorrectCatch' function because in
# for doing so the catch production function should be used.
#-------------------------------------------------------------------------------

CorrectCatch <- function(fleets, biols, BDs, biols.ctrl,fleets.ctrl, year = 1, season = 1,...){

    fleets <- unclass(fleets)
    yr <- year
    ss <- season
    it    <- dim(biols[[1]]@n)[6]
    nst <- length(biols)
    
    stnms <- names(biols)
    flnms <- names(fleets)             
    
    cth <-  matrix(fleets.ctrl$catch.threshold[,year,,season,drop=T],nst,it, 
                dimnames = list(dimnames(fleets.ctrl$catch.threshold)[[1]], 1:it)) # matrix[nstk,nit]

    Ba   <- lapply(stnms, function(x){   # biomass at age in the middle  of the season, list elements: [na,it] if age structured, [1,it] if biomass.
      if(dim(biols[[x]]@n)[1] > 1)
        return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss])
      else{
        if(biols.ctrl[[x]] == 'fixedPopulation'){
          return((biols[[x]]@n)[,yr,,ss, drop=F])
        }
        else{
          return((biols[[x]]@n + BDs[[x]]@gB)[,yr,,ss, drop=F])
        } }})
    names(Ba) <- stnms
    
    B    <- matrix(t(sapply(stnms, function(x){   # biomass in the middle if age struc. of the season  [ns,it]
      if(dim(biols[[x]]@n)[1] > 1)
        return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss, drop=T])
      else{
        if(biols.ctrl[[x]][['growth.model']] == 'fixedPopulation'){
          return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss, drop=T])
        }
        else{
          return((biols[[x]]@n*biols[[x]]@wt + BDs[[x]]@gB)[,yr,,ss, drop=T])
        }
        
      }})) , nst,it, dimnames = list(stnms, 1:it))
    for(st in stnms){
    
   #    print(st)

        if(dim(Ba[[st]])[1] > 1){ # age structured

            Cat  <- catchStock(fleets, st)[,yr,,ss]
            
            # Convert the [age,unit] combination into a continuous age.
            Ba.  <- unit2age(Ba[[st]]) # [na*nu,1,1,it]
            Cat. <- unit2age(Cat)
            K.   <- array(1,dim = dim(Ba.))   # Catch multipliers

            # CORRECT Ca if Ca > Ba, the correction is common for all the fleets.
            for(i in 1:it){

                if(any((Ba[[st]][,,,,,i]*cth[st,i] - Cat[,,,,,i]) < 0)){

                    cat('Ba*cth < Ca, for some "a" in stock',st, ', and iteration ', i,  '\n')

                    a.minus         <- which(Ba.[,,,i]*cth[st,i] < Cat.[,,,i])
                    a.plus          <- which(Ba.[,,,i]*cth[st,i] >= Cat.[,,,i])
                    K.[a.minus,,,i] <- Ba.[a.minus,,,i]*cth[st,i]/Cat.[a.minus,,,i]

                    # The correction below would correspond with a compensation of the decrease in 'a.minus' ages.
                    # K.[a.plus,,,i]  <- (Ct[i] - sum(Ba.[a.minus,,,i]))/sum(Cat.[a.plus,,,i])
                }
            }
            K <- age2unit(K., Ba[[st]])  # [na,1,nu,1,1,it]
            K[1,,-(1:ss)] <- 1   # This recruits do not exists  yet.
        }
        else{ # biomass dynamic.
        #   browser()
            Ct  <- c(catchWStock(fleets, st)[,yr,,ss,drop = F])  #[it]
            Bst <- c(array(B[st,], dim = c(1,1,1,1,1,it)))     # [it]
            K   <- rep(1,it)  # Catch multipliers

            if(any((Bst*cth[st,] - Ct) < 0)){
                i.minus <-  which((Bst*cth[st,] - Ct) < 0)
                cat('B*cth < C, for  stock',st, ', and iteration(s) ', i.minus,  '\n')
                K[i.minus] <- Bst[i.minus]*cth[st,i.minus]/Ct[i.minus]

            }

            K <- FLQuant(K, dimnames = dimnames(biols[[st]]@n[,yr,,ss]))
        }

        if(all(K==1)) next
        # Correct the catch fleet by fleet.
        # does the fleet_metier catch the stock??  which of them catch the stock??
        flinfo  <- stock.fleetInfo(fleets)
        flmtpos <-  colnames(flinfo)[which(flinfo[st, ] == 1)]

        for(k in flmtpos){
            k. <- strsplit(k, "&&")[[1]]
            fl <- k.[1]
            mt <- k.[2]

            cobj <- fleets[[fl]][[mt]][[st]]
            
     #       print(c(cobj@landings.n[14,yr,,ss]))
     

            cobj@landings.n[,yr,,ss] <-  cobj@landings.n[,yr,,ss]*K
            cobj@discards.n[,yr,,ss] <-  cobj@discards.n[,yr,,ss]*K
            cobj@landings[,yr,,ss]   <-  apply((cobj@landings.n*cobj@landings.wt)[,yr,,ss], c(2,4,6), sum)
            cobj@discards[,yr,,ss]   <-  apply((cobj@discards.n*cobj@discards.wt)[,yr,,ss], c(2,4,6), sum)

     #        print(c(cobj@landings.n[14,yr,,ss]))
             
            mt_idx <- which(names(fleets[[fl]]@metiers)==mt)
            fleets[[fl]] <- fill_flcatches(fl=fleets[[fl]],cobj=cobj,st=st,mt_idx=mt_idx)
        }
        
        
    }
    
    fleets <- FLFleetsExt(fleets)
    return(fleets)
}




                            
#-------------------------------------------------------------------------------
# CobbDouglasCom.CAA(fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
# CobbDouglasComb.CAA <- function(fleets, biols, BDs, biols.ctrl, fleets.ctrl, advice, year = 1, season = 1, flnm = 1, stknm = 1,...){
#   
#   if(dim(biols[[stknm]]@n)[1] > 1) return(CobbDouglasAge.CAA(fleets, biols, BDs, biols.ctrl, fleets.ctrl, advice, year, season, flnm, stknm))
#   if(dim(biols[[stknm]]@n)[1] == 1) return(CobbDouglasBio.CAA(fleets, biols, BDs, biols.ctrl, fleets.ctrl, advice, year, season, flnm, stknm))
# }
# 


#-------------------------------------------------------------------------------
# Baranov (fleets, biols, year = 1, season = 1)
#-------------------------------------------------------------------------------
# Baranov.CAA <- function(fleets, biols, BDs, biols.ctrl, fleets.ctrl, advice, year = 1, season = 1, flnm = 1, stknm = 1,...){
#   
#   rho <- fleets.ctrl[['catch.threshold']][stknm,year,, season,drop=T] # [it]
#   
#   nf    <- length(fleets)
#   stnms <- names(biols)
#   nst   <- length(stnms)
#   it    <- dim(biols[[1]]@n)[6]
#   
#   #   if(year == 35 & stknm == 'HKE') browser()
#   
#   fleets <- unclass(fleets)
#   
#   yr <- year
#   ss <- season
#   st <- stknm
#   
#   fl    <- fleets[[flnm]]
#   sts   <- catchNames(fl)
#   mtnms <- names(fl@metiers)
#   
#   if(!(st %in% sts)) return(fleets)
#   
#   # tac <- rep(Inf,it)
#   
#   # catch restriction, if empty => landings.
#   if (is.null(fleets.ctrl[[flnm]]$restriction)) {
#     catch.restr <- 'landings'
#   } else 
#     catch.restr <- ifelse(length(fleets.ctrl[[flnm]]$restriction)==1, fleets.ctrl[[flnm]]$restriction, 
#                           fleets.ctrl[[flnm]]$restriction[yr])
#   # catch.restr <- ifelse(is.null(fleets.ctrl[[flnm]]$restriction), 'landings',ifelse(length(fleets.ctrl[[flnm]]$restriction)==1, fleets.ctrl[[flnm]]$restriction,fleets.ctrl[[flnm]]$restriction[yr]))
#   
#   
#   #  quota share % to be upodate du to year transfer.
#   yrtr_p <- fleets.ctrl[[flnm]]$LandObl_yearTransfer_p[stknm,yr]
#   yrtr_p <- ifelse(is.null(yrtr_p), 0,yrtr_p)
#   # if year transfer was used in previous year discount it, absolute catch
#   yrtr_disc <- fleets.ctrl[[flnm]]$LandObl_discount_yrtransfer[stknm,yr-1,] # [it]
#   
#   if(is.null(yrtr_disc)) yrtr_disc <- 0
#   
#   fleets.ctrl[[flnm]] 
#   # if TAC overshoot is discarded, calculate seasonal TAC to calculate the discards.
#   TACOS <- fleets.ctrl[[flnm]][[stknm]][['discard.TAC.OS']]    # Is the TAC overshot discarded?
#   TACOS <- ifelse(is.null(TACOS), TRUE, TACOS) 
#   
#   yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
#   ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
#   QS          <- yr.share*ss.share                                          # [it]
#   QS[is.na(QS)] <- 0              
#   tac <- (advice$TAC[st,yr,drop=T])*QS 
#   
#   #in case there is year transfer the tac is recalculated
#   
#   if(TACOS){
#     yr.share    <- advice$quota.share[[stknm]][flnm,yr,, drop=T]              # [it]
#     ss.share    <- fleets.ctrl$seasonal.share[[stknm]][flnm,yr,,ss, drop=T]   # [it]
#     QS          <- yr.share*ss.share                                          # [it]
#     QS[is.na(QS)] <- 0              
#     tac <- (advice$TAC[st,yr,drop=T]*QS*(1+yrtr_p)) - yrtr_disc # it, add yeartransfer in case it is in place, first we increment in % the quota and then we discount the cuota used in previous year. 
#     # the minimise is not added because it is discarded.      
#   }
#   
#   
#   # Re-scale QS to fleet share within the season instead of season-fleet share within year.
#   # list of stocks, each stock [nf,it]
#   
#   
#   if(dim(biols[[st]]@n)[1] == 1) stop(st, ' stock has no ages, Cobb Douglas cannot be applied at age level then! correct the "catch.model" argument in "fleets.ctrl" argument!\n')
#   
#   N <- (biols[[st]]@n*exp(-biols[[st]]@m/2))[,yr,,ss]  # Ba[na,it], biomass at age in the middle  of the season,
#   B <- unitSums(quantSums(N*biols[[st]]@wt[,yr,,ss]))
#   
#   #TAC.SS is total catch in the season that is going to used to estimate total F, in CBaranov,
#   #for that we need to estimate the QS.SS
#   
#   yr.share.All    <- matrix(advice$quota.share[[stknm]][,yr,, drop=T],
#                             nrow=dim(advice$quota.share[[stknm]][,yr,, ])[1],ncol=it)        # [nf,it]
#   ss.share.All    <- matrix( fleets.ctrl$seasonal.share[[stknm]][,yr,,ss, drop=T]  ,
#                              nrow=dim(advice$quota.share[[stknm]][,yr,, ])[1],ncol=it)  
#   # [nf,it]
#   QS.ss <- apply(yr.share.All*ss.share,2,sum ) #[it]
#   QS.ss[is.na(QS.ss)] <- 0
#   
#   
#   TAC.yr  <- advice$TAC[stknm,yr,drop=T]
#   
#   #  for(stknm in  stnms){
#   tacos.fun <- fleets.ctrl[[flnm]][[stknm]]$TAC.OS.model
#   if(is.null(tacos.fun))   alpha <- rep(1,it)
#   else{
#     alpha <- eval(call(tacos.fun, fleets = fleets, TAC = TAC.yr, fleets.ctrl = fleets.ctrl, flnm = flnm, stknm = stknm, year = year, season = season))
#   }
#   TAC.yr <- TAC.yr*alpha 
#   # }
#   
#   TAC.ss<- TAC.yr*QS.ss
#   TAC.ss<- ifelse(as.vector(B*rho) < TAC.ss, as.vector(B*rho), TAC.ss)
#   Cr.f <- tac #catches per fleet
#   
#   efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])),
#                   length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
#   eff   <- matrix(fl@effort[,yr,,ss],length(mtnms), it, dimnames = list(mtnms, 1:it), byrow = T)
#   
#   # flinfo: matrix with information on which metier catch which stock.
#   fl.        <- FLFleetsExt(fl)
#   names(fl.) <- flnm
#   flinfo     <- stock.fleetInfo(fl.)
#   flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')
#   
#   mtst <- flinfo[[st]][2]
#   
#   age.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[1]]
#   age.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[1]]
#   age.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[1]]
#   
#   unit.q     <- dimnames(fl@metiers[[mtst]]@catches[[st]]@catch.q)[[3]]
#   unit.alpha <- dimnames(fl@metiers[[mtst]]@catches[[st]]@alpha)[[3]]
#   unit.beta  <- dimnames(fl@metiers[[mtst]]@catches[[st]]@beta)[[3]]
#   
#   q.m   <- array(0, dim = c(length(mtnms), length(age.q), length(unit.q),it),     dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1:it))
#   alpha.m <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), it), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1:it))
#   beta.m  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
#   ret.m  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
#   wl.m   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
#   wd.m   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
#   
#   
#   for(mt in mtnms){
#     
#     if(!(st %in% names(fl@metiers[[mt]]@catches))) next
#     
#     q.m[mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE]
#     alpha.m[mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE]
#     beta.m[mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE]
#     ret.m[mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, drop = TRUE]
#     wl.m[mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, drop = TRUE]
#     wd.m[mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, drop = TRUE]
#   }
#   
#   Nst  <- array(N[drop=T],dim = dim(N)[c(1,3,6)])
#   Cam <- array(NA, dim=dim(q.m))
#   
#   for(i in 1: it){
#     
#     Cayr_1 <- catchStock(fleets,st)[,yr-1,,,,i,drop=F]
#     Nayr_1 <- biols[[st]]@n[,yr-1,,ss,,i,drop=F]
#     Mayr_1 <- biols[[st]]@m[,yr-1,,ss,,i,drop=F]
#     Na  <-  biols[[st]]@n[,yr,,ss,,i,drop=F] # N is in the middle of the season,
#     #and we need N at athe beginning of the season
#     Ma <- biols[[st]]@m[,yr,,ss,,i,drop=F]
#     Wa <- biols[[st]]@wt[,yr,,ss,,i,drop=F]
#     n.mt <- length(mtnms)
#     Cafyr_1 <- Nayr_1
#     Cafyr_1[] <- 0
#     
#     for (met in 1:n.mt){
#       Cafyr_1 <- Cafyr_1+iter((fleets[[flnm]]@metiers[[met]]@catches[[st]]@landings.n+
#                                  fleets[[flnm]]@metiers[[met]]@catches[[st]]@discards.n)[,yr-1,,ss],it)}
#     
#     Cam[,,,i] <- Baranov(E = eff[1,], fleets=fleets,biols=biols, Cr=Cr.f[i],N = Nst[,,i],   efs.m = efs.m[,i,drop=FALSE], q.m = q.m[,,,i,drop=FALSE], 
#                           alpha.m = alpha.m[,,,i], beta.m = beta.m[,,,i], wd.m =wd.m [,,,i,drop=FALSE],
#                           wl.m = wl.m[,,,i,drop=FALSE], ret.m = ret.m[,,,i,drop=FALSE],
#                           tac=TAC.ss[i],Cayr_1=Cayr_1,Nayr_1=Nayr_1,Mayr_1=Mayr_1,Na=Na,Ma=Ma,Cafyr_1=Cafyr_1)
#   }
#   # Cam <- CobbDouglasAge(E = eff[1,], N = Nst, wl.m = wl.m, wd.m = wd.m, ret.m = ret.m, q.m = q.m,
#   #                        efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, rho = rho)
#   
#   
#   #    if(st == 'LDB' & flnm == 'DTS_SP') browser()
#   
#   # if catch restriction is landings, Lrat is calculated over landigns, else it is calculated over total catch including undersize individuals.
#   Ctotal <- ifelse(rep(catch.restr == 'landings', it), apply(Cam*ret.m,4,sum), apply(Cam,4,sum)) 
#   
#   tac.disc <- ifelse(Ctotal < tac, rep(1,it), tac/Ctotal)
#   tac.disc <- ifelse(Ctotal == 0, 0, tac.disc) # In case Ctotal= 0 to avoid NaN in tac.disc
#   
#   
#   # cat('Lrat: ', tac.disc, '\n')
#   # cat('C: ', Ctotal, '\n')
#   
#   Cam <- array(Cam, dim = c(length(mtnms),dim(biols[[st]]@n)[1], 1, dim(biols[[st]]@n)[3],1,it))
#   
#   for(mt in 1:length(mtnms)){
#     
#     Ca <- array(Cam[mt,,,,,], dim = c(dim(biols[[st]]@n)[1], 1, dim(biols[[st]]@n)[3],1,1,it))
#     
#     if(!(st %in% names(fl@metiers[[mt]]@catches))) next
#     cobj <- fl[[mt]][[st]]
#     
#     
#     na <- dim(q.m)[2]
#     nu <- ifelse(is.na(dim(q.m)[3]), 1, dim(q.m)[3])
#     
#     efm <- array(eff[1,]*efs.m[mt,], dim = c(it,na,1,nu,1,1))
#     efm <- aperm(efm, c(2:6,1))
#     
#     dsa <- cobj@discards.sel[,yr,,ss]  # [na,1,nu,1,1,it]
#     lsa <- cobj@landings.sel[,yr,,ss]  # [na,1,nu,1,1,it]
#     sa  <- (dsa + lsa)  
#     
#     # Recalculate dsa and lsa according to 'tac.disc'     # [na,nu,it]
#     lsa <- sweep(lsa,6,tac.disc, "*")
#     dsa <- sa - lsa             # [na,nu,it]
#     
#     cobj@discards.n[,yr,,ss] <- Ca*dsa/sa/cobj@discards.wt[,yr,,ss]
#     cobj@landings.n[,yr,,ss] <- Ca*lsa/sa/cobj@landings.wt[,yr,,ss]
#     
#     # When sa = 0 <-  land.n & dis.n = NA => change to 0.
#     cobj@landings.n[,yr,,ss][is.na(cobj@landings.n[,yr,,ss])] <- 0
#     cobj@discards.n[,yr,,ss][is.na(cobj@discards.n[,yr,,ss])] <- 0
#     
#     cobj@discards[,yr,,ss] <- apply(Ca*dsa/sa,c(2,4,6),sum,na.rm=T)
#     cobj@landings[,yr,,ss] <- apply(Ca*lsa/sa,c(2,4,6),sum,na.rm=T)
#     
#     fl@metiers[[mt]]@catches[[st]] <- cobj          
#   }
#   
#   fleets[[flnm]] <- fl
#   
#   #    fleets <- FLFleetsExt(fleets)
#   
#   return(fleets)
# }
