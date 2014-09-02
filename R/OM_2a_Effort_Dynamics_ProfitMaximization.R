#-------------------------------------------------------------------------------
# Functions to be used within Rdonlp2 to obtain effort and effortshare according
# to maximization of profits, restringed to Fcube like approach in 
# TAC constraints.
#
# Dorleta Garcia - Azti Tecnalia
# Created: 22/11/2011 19:02:56
# Changed: 23/11/2011 10:16:01
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# MaxProfit.stkCnst(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1)
# Catch production function based on Cobb-Doug at age level.
# OVER-QUOTA LANDINGS ARE DISCARDED => NOT BENEFIT FROM THEM.
#-------------------------------------------------------------------------------
MaxProfit.stkCnst <- function(fleets, biols, covars, advice, fleets.ctrl, flnm, year = 1, season = 1){
    
    dimnms <- dimnames(biols[[1]]@n)
    
    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )

    # If year/season/iter numerics => indicate position 
    # else names => get positions.
    if(length(year) > 1 | length(season) > 1)
        stop('Only one year and season is allowed' )
    
    # 'year' dimension.
    yr <- year
    if(is.character(year)) yr <- which(dimnms[[2]] %in% year)
    if(length(yr) == 0) stop('The year is outside object time range')  
    
    # 'season' dimension.
    ss <- season
    if(is.character(season)) ss <- which(dimnms[[4]] %in% season)
    if(length(ss) == 0) stop('The season is outside object season range')  
    
    # Check fleets.ctrl elements.
 # This checks the restriction for all the fleets and we only need to check the fleet we are proyecting.
 #   if(! all(sapply(names(fleets), function(x) fleets.ctrl[[x]]$restriction %in% c('catch', 'landings'))))
 #       stop("fleets.ctrl$restriction must be equal to 'catch' or 'landings'")
    if(!(fleets.ctrl[[flnm]]$restriction %in% c('catch', 'landings') )) stop("fleets.ctrl$restriction for fleet, ', flnm, ', must be equal to 'catch' or 'landings'")

    # Dimensions.
    nst <- length(biols);          stnms <- names(biols)
    ns  <- dim(biols[[1]]@n)[4]
    it  <- dim(biols[[1]]@n)[6]
    flnms <- names(fleets)
    

        
    # Biomass at age.
    B    <- matrix(t(sapply(stnms, function(x){   # biomass in the middle of the season  [nst,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return(unitSums(quantSums(biols[[x]]@n*biols[[x]]@wt*exp(-biols[[x]]@m/2)))[,yr,,ss, drop=T])
                                else return((biols[[x]]@n*biols[[x]]@wt)[,yr,,ss, drop=T])})) , nst,it, dimnames = list(stnms, 1:it))

    N   <- lapply(stnms, function(x){   # biomass at age in the middle  of the season, list elements: [na,1,nu,1,1,it]
                                if(dim(biols[[x]]@n)[1] > 1)
                                    return((biols[[x]]@n*exp(-biols[[x]]@m/2))[,yr,,ss, drop = FALSE])
                                else return((biols[[x]]@n)[,yr,,ss])})
    names(N) <- stnms

                          
    QS.fls   <- lapply(stnms, function(x){           # list of stocks, each stock [nf,it]
                            # Calculate QS by fleet for the year and season
                            yr.share    <- advice$quota.share[[x]][,yr,, drop=T]        # [nf,it]
                            ss.share    <- fleets.ctrl$seasonal.share[[x]][,yr,,ss, drop=T]   # [nf,it]
                            quota.share <-  matrix(yr.share*ss.share, length(flnms), it, dimnames = list(flnms, 1:it))
                            quota.share[is.na(quota.share)] <- 0
                            return(quota.share)})         
    names(QS.fls) <- stnms
    
        
    # If TAC >= B*alpha => TAC = B*alpha.
    TAC.yr   <- matrix(advice$TAC[,yr,drop=T], nst,it, dimnames = list(stnms, 1:it))                       # [nst,it]
    CT       <- fleets.ctrl$catch.threshold[,yr,,ss, drop=T]  # [ns,it]
    QS.ss    <- matrix(t(sapply(stnms, function(x) apply(QS.fls[[x]],2,sum))), nst,it, dimnames = list(stnms, 1:it))  # [nst,it]
                            
    TAC <- ifelse(B*CT < TAC.yr*QS.ss, B*CT, TAC.yr*QS.ss) 
    
    # Re-scale QS to fleet share within the season instead of season-fleet share within year.
    QS   <- lapply(stnms, function(x){          # list of stocks, each stock [nf,it]
                            res <- sweep(QS.fls[[x]], 2, apply(QS.fls[[x]],2, sum), "/")
                            res[is.na(res)] <- 0 
                            return(res)})      
    names(QS) <- stnms

    fl    <- fleets[[flnm]]
    sts   <- catchNames(fl)
    mtnms <- names(fl@metiers)
    nmt   <- length(mtnms)
    
    # flinfo: matrix with information on which metier catch which stock.
    fl.        <- FLFleetsExt(fl)
    names(fl.) <- flnm
    flinfo     <- stock.fleetInfo(fl.)
    flinfo <-  strsplit(apply(flinfo, 1,function(x) names(which(x == 1))[1]), '&&')

    efs.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@effshare[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    vc.m <- matrix(t(sapply(mtnms, function(x) fl@metiers[[x]]@vcost[,yr,,ss, drop=T])), 
                    length(mtnms), it, dimnames = list(metier = mtnms, 1:it))
    fc    <- fl@fcost[,yr,,ss, drop=T]* covars$NumbVessels[flnm,yr,,ss, drop=T] # [it]
    crewS <- fl@crewshare[,yr,,ss, drop=T] # [it]
    if(length(fc) == 1) fc <- rep(fc,it)
    
    effs <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))
    Cr.f <- matrix(NA,length(sts), it, dimnames = list(sts, 1:it))
        

    q.m <- alpha.m <- beta.m  <- pr.m <- ret.m <- wd.m <- wl.m <-vector('list', length(sts))
    names(q.m) <- names(pr.m) <- names(alpha.m) <- names(beta.m) <- names(ret.m) <- names(wl.m) <- names(wd.m) <- sts

    tacos <- logical(length(sts))
    names(tacos) <- sts
    
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
            
        q.m[[st]]     <- array(0, dim = c(length(mtnms), length(age.q),     length(unit.q),it),      dimnames = list(metier = mtnms, age = age.q, unit = unit.q, iter = 1:it))
        alpha.m[[st]] <- array(0, dim = c(length(mtnms), length(age.alpha), length(unit.alpha), it), dimnames = list(metier = mtnms, age = age.q, unit = unit.alpha, iter = 1:it))
        beta.m[[st]]  <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        ret.m[[st]]   <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        wl.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        wd.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))
        pr.m[[st]]    <- array(0, dim = c(length(mtnms), length(age.beta), length(unit.beta), it),  dimnames = list(metier = mtnms, age = age.beta,unit = unit.beta,  iter = 1:it))


        # if TAC overshoot is not discarded, it is sold and it contributes to the revenue.
        tacos[st] <- ifelse(is.null(fleets.ctrl[[flnm]][[st]][['discard.TAC.OS']]), TRUE,fleets.ctrl[[flnm]][[st]][['discard.TAC.OS']]) 
        
        for(mt in mtnms){
            if(!(st %in% names(fl@metiers[[mt]]@catches))) next
                    
            q.m[[st]][mt,,,]     <- fl@metiers[[mt]]@catches[[st]]@catch.q[,yr,,ss, drop = TRUE] 
            alpha.m[[st]][mt,,,] <- fl@metiers[[mt]]@catches[[st]]@alpha[,yr,,ss, drop = TRUE] 
            beta.m[[st]][mt,,,]  <- fl@metiers[[mt]]@catches[[st]]@beta[,yr,,ss, drop = TRUE] 
            ret.m[[st]][mt,,,]   <- fl@metiers[[mt]]@catches[[st]]@landings.sel[,yr,,ss, drop = TRUE] 
            wl.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@landings.wt[,yr,,ss, drop = TRUE]
            wd.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@discards.wt[,yr,,ss, drop = TRUE]

            # The price is taken from the year before, because price for the current year is updated after catch is produced,
            # if the price was dynamically updated inside this function the optimizer could crash. 
            pr.m[[st]][mt,,,]    <- fl@metiers[[mt]]@catches[[st]]@price[,yr-1,,ss, drop = TRUE] 
        }    
        Cr.f[st,] <- TAC[st,]*QS[[st]][flnm,]
    }

    Et  <- fl@effort[,yr-1,drop=T]   # it
    efs <- t(matrix(sapply(fl@metiers, function(x) x@effshare[,yr,,ss,drop=T]),it,nmt)) #[nmt, it]   
    
    Et.res  <- numeric(it)
    efs.res <- matrix(NA,nmt,it, dimnames = list(mtnms,1:it))
        
    # CALCULATE THE EFFORT THAT MAXIMIZES THE PROFITS UNCER CERTAIN CONSTRAINTS.
    # (Rdonlp2 library)    
    stk.cnst <- fleets.ctrl[[flnm]]$stk.cnst
    K <- fl@capacity[,yr,,ss,drop=T]
    if(length(K) == 1) K <- rep(K,it)
        
    ctrl <- donlp2Control() #,
    ctrl$fnscale <-  -1 
    ctrl$silent  <- TRUE
    #ctrl$report <-  TRUE 
    #ctrl$epsx  <-  1e-5
    #ctrl$epsdif <- 1e-4      
        
    res <- vector('list', it)
    for(i in 1:it){
        
        donlpDat <- list(q.m = q.m, alpha.m = alpha.m , beta.m = beta.m, pr.m = pr.m, ret.m = ret.m, wl.m = wl.m, wd.m = wd.m, crewS = crewS,
                         vc.m = vc.m, N = N, fc = fc, stk.cnst = stk.cnst, Cr.f = Cr.f, i = i,
                         catch.rest = fleets.ctrl[[flnm]]$restriction, tacos = tacos)
        attach(donlpDat)
            
        E <- rep(K[i]/nmt,nmt)*0.5 # c(Et[i]*efs[,i])   
        #      save(donlpDat,E,K, file = 'inputData.Rdata')
        #       browser()
        res[[i]] <- donlp2(E, fobj.maxprofits, par.upper = rep(Inf,nmt), par.lower = rep(.1, nmt),  # Obj. Function and params. bound
                     A = matrix(rep(1,nmt),1,nmt), lin.upper = K[i], lin.lower = .1*nmt, control = ctrl,#,                      # Linear constraints
                     nlin = list(nlin.maxprofits), nlin.upper =  Cr.f[stk.cnst,i], nlin.lower = 0)         # Non-Linear constraints
        Et.res[i]   <- sum(res[[i]]$par)
        efs.res[,i] <- res[[i]]$par/sum(res[[i]]$par)
        cat('Effort share: ', efs.res[,i], ', ~~~~~ Effort: ',Et.res[i], ', ~~~~~ Benefit: ', res[[i]]$fx, '\n')  
            detach()
    }
        
        
    fleets[[flnm]]@effort[,yr,,ss] <- Et.res
    for(mt in mtnms)  fleets[[flnm]]@metiers[[mt]]@effshare[,yr,,ss] <- efs.res[mt,]        
                    
    # Update the quota share of this step and the next one if the 
    # quota share does not coincide with the actual catch. (update next one only if s < ns).
    for(st in sts){
      #   browser()
      #  if(st == 'SKH') browser()
        
        yr.share       <- advice$quota.share[[st]][flnm,yr,, drop=T]      # [it]
        ss.share       <- t(matrix(fleets.ctrl$seasonal.share[[st]][flnm,yr,,, drop=T], ns, it))# [it,ns]
        quota.share.OR <- matrix(t(yr.share*ss.share), ns, it)
        # The catch.
        catchFun <- paste(fleets.ctrl[[flnm]][[st]][['catch.model']], 'CatchFleet', sep = ".")
        Dim <- dim(N[[st]])
        N. <- array(dim = c(Dim[1],1,Dim[2],1,1,Dim[6]))
        N.[,1,,1,1,]  <- N[[st]]
        catch <- eval(call(catchFun, N = N., B = B[st,], effort = Et.res, efs.m = efs.res, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]], wl.m = wl.m[[st]], wd.m = wd.m[[st]], ret.m = ret.m[[st]]))
            
        quota.share    <- updateQS.SMFB(QS = quota.share.OR, TAC = TAC.yr[st,], catch = catch, season = ss)        # [ns,it]
                              
        fleets.ctrl$seasonal.share[[st]][flnm,yr,,] <- t(t(quota.share)/apply(quota.share, 2,sum)) #[ns,it], doble 't' to perform correctly de division between matrix and vector.
    }

    
    return(list(fleets = fleets, fleets.ctrl = fleets.ctrl, res = res))
}



#-------------------------------------------------------------------------------
# fobj.maxprofits(E) :: Objective function 
#                       (function to be maximized in a fleet by fleet case)
#    - E: numeric(nmt + 1). E[1]  = Total effort
#                           E[-1] = Effort shares across metiers.
#    - qa.mt ~ alpha.a.mt ~ beta.a.mt ~ pra.mt.i :: lists with one element per
#                   stock, each element with the form: array(na_st,nmt,it)
#    - Ba ~ list with one element per stock, each element with the form: array(na_st,it)
#                       
# NOTE: The objects qa.mt, alpha.a.mt, beta.a.mt, pra.mt.i & Ba _MUST_ be in the 
#    global environment when this function is used in "rdonlp2" function. The 
#    reason is that the function does not allow extra arguments for the 
#    non-control variables.
#    Furthermore, 'i' (the iteration) must also be in the global environment, 
#    in order to optimize the computations the data for all the iterations is
#    extract at once, an then each iteration's data is extracted within the 
#    function.     
#-------------------------------------------------------------------------------
        
fobj.maxprofits <- function(E){


    nmt <- length(E)

    # Data, extract the iterations.
    q.m.i      <- lapply(q.m, function(x) x[,,,i, drop=F])      # [nmt,na,nu,1]
    alpha.m.i  <- lapply(alpha.m, function(x) x[,,,i, drop=F])  # [nmt,na,nu,1]
    beta.m.i   <- lapply(beta.m, function(x) x[,,,i, drop=F])   # [nmt,na,nu,1]
    pr.m.i     <- lapply(pr.m, function(x) x[,,,i, drop=F])     # [nmt,na,nu,1]
    ret.m.i     <- lapply(ret.m, function(x) x[,,,i, drop=F])     # [nmt,na,nu,1]
    wd.m.i     <- lapply(wd.m, function(x) x[,,,i, drop=F])     # [nmt,na,nu,1]
    wl.m.i     <- lapply(wl.m, function(x) x[,,,i, drop=F])     # [nmt,na,nu,1]
    vc.m.i     <- vc.m[,i] # [nmt]
    N.i       <- lapply(names(q.m.i), function(x){          # [nmt,na,nu]
                                    Dim <- dim(q.m.i[[x]])
                                    res <- array(dim = Dim)
                                    for(j in 1:nmt) res[j,,,] <- N[[x]][,,,,,i,drop=F]
                                    return(res)})
    Cr.f.i <- Cr.f[,i,drop=F] # [nst,1]   quota share 
    
    res <- 0
    

    # Cst <- rep(0,4)
    # names(Cst) <- c('NHKE', 'CMEG', 'CANK', 'CMON')
    
    Cst <- Lst <-  NULL
    
    for(st in 1:length(q.m.i)){
    

        E  <- array(E,dim = dim(N.i[[st]]))   # [nmt,na,nu]

        Cst[st] <- sum(q.m.i[[st]]*(E^alpha.m.i[[st]])*
                (N.i[[st]]*(ret.m.i[[st]]*wl.m.i[[st]]+ (1-ret.m.i[[st]])*wd.m.i[[st]])^beta.m.i[[st]]))
        Lst[st] <- sum(ret.m.i[[st]]*q.m.i[[st]]*(E^alpha.m.i[[st]])*
                (N.i[[st]]*(ret.m.i[[st]]*wl.m.i[[st]]+ (1-ret.m.i[[st]])*wd.m.i[[st]])^beta.m.i[[st]]))  # multiply the retention vector if landing is the restriction.
        
        # The oversized discards are always discarded, but if landing obligation is in place they account in quota (catch == TRUE). 
        if(catch.rest  != 'catch') Cst[st] <- Lst[st]# The restriction is landings.
        
        # TAC overshot can be landed or discarded. In the case of landing obligation it is 'discarded' because it does not 
        # contribute to the revenue but it goes against the TAC => TACOS == TRUE
        if(tacos[st] == FALSE) # TAC Overshot is not discarded.
            Lrat <- 1
        else Lrat <- ifelse(Cr.f.i[st,]/Cst[st] > 1, 1, Cr.f.i[st,]/Cst[st])  # TAC Overshot is  discarded.
                                                                              # The overquota discards are proportional to the catch in all the metiers.
        
        
        res <- res + sum(q.m.i[[st]]*(E^alpha.m.i[[st]])*
                (N.i[[st]]*(ret.m.i[[st]]*wl.m.i[[st]]+ (1-ret.m.i[[st]])*wd.m.i[[st]])^beta.m.i[[st]]))*Lrat
     #   print(res)

    }

    res <- res - sum(vc.m.i*E[,1,1,1]) - fc[i]  - res*crewS[i]
#    print(res)
#    cat('--- NHKE: ', Cst[1], ' --- CMEG: ', Cst[2], ' --- CANK: ', Cst[3], ' --- CMON: ', Cst[4], '\n') 
#    cat('****** ', E[1:6], ' *******\n')
    return(res)
}

#-------------------------------------------------------------------------------
# nlin.maxprofits(E) 
#   * The maximization of profits is restringed to the compliance, 
#     in some degree, of the TACs. There will be different options based on the 
#     Fcube like approach. The contraint for all the stock will have similar form
#     so we write the general function  'nlin.maxprofits' and then we will use
#     it in accordance with the option choosen.
#-------------------------------------------------------------------------------
nlin.maxprofits <- function(E){

    nmt <- length(E)
    na  <- dim(N[[stk.cnst]])[1]

    # Data, extract the iterations.
    q.m.i      <- q.m[[stk.cnst]][,,,i, drop=F]      # [nmt,na,nu,1]
    alpha.m.i  <- alpha.m[[stk.cnst]][,,,i, drop=F]  # [nmt,na,nu,1]
    beta.m.i   <- beta.m[[stk.cnst]][,,,i, drop=F]   # [nmt,na,nu,1]
    wd.m.i     <- wd.m[[stk.cnst]][,,,i, drop=F]   # [nmt,na,nu,1]
    wl.m.i     <- wl.m[[stk.cnst]][,,,i, drop=F]   # [nmt,na,nu,1]
    ret.m.i    <- ret.m[[stk.cnst]][,,,i, drop=F]   # [nmt,na,nu,1]
    pr.m.i     <- pr.m[[stk.cnst]][,,,i, drop=F]     # [nmt,na,nu,1]
    N.i       <- E1 <- array(dim = dim(q.m.i))

    for(j in 1:nmt) N.i[j,,,] <- N[[stk.cnst]][,,,,,i,drop=F]   # [nmt,na,nu,1]

    E1  <- array(E, dim(q.m.i))

    res <-    sum(q.m.i*(E^alpha.m.i)*(N.i*(ret.m.i*wl.m.i + (1-ret.m.i)*wd.m.i)^beta.m.i))

    return(res)
}


# MIN OPTION
#   * The maximization of benefits is constrained to the TACs of all the stocks
#     (if one stock is not managed put each TAC share equal to infinity). 
#     That is, the catch must be below TAC quota share.
#   * In order to do that we 
#-------------------------------------------------------------------------------




