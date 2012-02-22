#------------------------------------------------------------------------------#
#              Tools to help obtaining an initial random population
#
# XSAboot - Generates a reandom population bootstraping XSA's indices' residuals.
#
# 08/02/2011 14:26:20
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#  XSAboot(Stock, Indices, control, parametric = TRUE)
# - Function taken from FLMSSim used extensively in NHake LTMP (2007).
#------------------------------------------------------------------------------#

XSAboot <- function(Stock, Indices, control, parametric = TRUE, mvnorm = TRUE){

    it <- dim(Stock@m)[6]
    na <- dim(Stock@m)[1]
    ny <- dim(Stock@m)[2]
    XSA <- vector('list', it)

    # Run an XSA with the first iteration of Stock and Indices.

    stk  <- qapply(Stock,E)
    inds <- FLIndices(lapply(Indices,function(x) qapply(x,E)))
    
    xsa.base <- FLXSA(stk, inds, control, diag.flag = TRUE)

    ind    <- xsa.base@index
    logres <- xsa.base@index.res
    q.hat  <- xsa.base@q.hat

    dms <- lapply(logres, dim)
    dmnms <- lapply(logres, dimnames)
    varresids <- lapply(xsa.base@index.res, function(x) c(apply(x, 1, var, na.rm = TRUE)))

    Z <- xsa.base@harvest + Stock@m

    if(parametric == TRUE){
        for(i in 1:length(ind)){

            Z. <- Z[dmnms[[i]]$age, dmnms[[i]]$year]
            alpha <- Indices[[i]]@range[['startf']]
            beta <- Indices[[i]]@range[['endf']]
            tau <- (exp(-alpha*Z.) - exp(-beta*Z.))/(Z.*(beta-alpha))

            # residuals in log-scale
            logres     <- xsa.base@index.res[[i]][drop = TRUE]
            # correlation at age of log-residuals
            cor_logres <- cor(t(xsa.base@index.res[[i]][drop = TRUE]))
            # stdev at age of log-residuals
            std_logres <- sqrt(apply(logres,1, var))
            # covariance at age of log-residuals
            cov_logres <- cor_logres
            na <- dim(cor_logres)[1]
            for(j in 1:na) for(k in 1:na) cov_logres[j,k] <- cor_logres[j,k]*(prod(std_logres[c(j,k)]))
            
            if(mvnorm == FALSE){
              for(j in 1:na) cov_logres[-j,-j] <- 0           
            }
            
            res0 <-  t(rmvnorm(it*ny,rep(0,na),cov_logres, method = 'svd'))   # [na,it*ny]
            res  <-  Indices[[i]]@index


            Indices[[i]]@index.var[] <- varresids[[i]]
            res <- apply(sqrt(Indices[[i]]@index.var), 1:6, rnorm, n = 1, mean = 0)
            Indices[[i]]@index.q <- FLQuant(c(q.hat[[i]]), dimnames = list(age = dmnms[[i]]$age, year = dmnms[[i]]$year), iter = it, fill.iter = TRUE)
            Indices[[i]]@index[] <- Indices[[i]]@index.q*xsa.base@stock.n[dmnms[[i]]$age, dmnms[[i]]$year]*tau*exp(res)
        }
    }
    else{ # non-parametric bootstrap of the residuals.
      logresids <- list()
      
        for(i in 1:length(ind)){
            logresids[[i]] <- Indices[[i]]@index
            Z. <- Z[dmnms[[i]]$age, dmnms[[i]]$year]
            alpha <- Indices[[i]]@range[['startf']]
            beta <- Indices[[i]]@range[['endf']]
            tau <- (exp(-alpha*Z.) - exp(-beta*Z.))/(Z.*(beta-alpha))

            for(a in 1:length(dimnames(Indices[[i]]@index)[[1]])){
                logresids[[i]][a,,] <- sample(c(logres[[i]][a,]), it*dms[[i]][2], replace = TRUE)}

             q.hat <- FLQuant(c(q.hat[[i]]), dimnames = list(age = dmnms[[i]]$age, year = dmnms[[i]]$year), iter = it, fill.iter = TRUE)
             Indices[[i]]@index[] <- q.hat*xsa.base@stock.n[dmnms[[i]]$age, dmnms[[i]]$year]*tau*exp(logresids[[i]])
        }
    }

     res  <- FLXSA(Stock, Indices, control, diag.flag = FALSE)

     Stock <- Stock + res
    

       return(list(Stock = Stock, Indices = Indices))
}


#------------------------------------------------------------------------------#
#  XSAboot(Stock, Indices, control, parametric = TRUE)
# - Function taken from FLMSSim used extensively in NHake LTMP (2007).
#------------------------------------------------------------------------------#
# Generate Stock
# Based upon the biological characteristics and the selection pattern it is possible to describe the expected equilibrium dynamics

initPop_lifeHist <-function(k,Linf,steepness=0.75,vbiomass=1e3,srModel="bevholt",sel=c(a=1,sL=1,sR=1),mat95=3,age=1:75,fmsy=rep(1.0,101))
   {
   ## Dimension stuff
   yrs<-1:length(fmsy)
   dms<-list(age=age,year=1)


   ## Biological stuff
   wts       <-FLQuant(vbMass(Linf,k,age),dimnames=dms)
   m         <-FLQuant(M(c(1000*vbMass(Linf,k,age+0.5))^(1/3),Linf),dimnames=dms)
   mat50     <-Mat50(c(mean(m)),k)
   selPattern<-dnormal(age,(mat50+mat95)*sel["a"],(mat50+mat95)*sel["sL"],(mat50+mat95)*sel["sR"])

   ## create FLBRP object
   res<-FLBRP(stock.wt       =wts,
              landings.wt    =wts,
              discards.wt    =wts,
              bycatch.wt     =wts,
              m              =m,
              mat            =FLQuant(logistic(age,mat50,mat95), dimnames=dms),
              landings.sel   =FLQuant(selPattern,dimnames=dms),
              discards.sel   =FLQuant(0,dimnames=dms),
              bycatch.harvest=FLQuant(0,dimnames=dms),
              harvest.spwn   =FLQuant(0,dimnames=dms),
              m.spwn         =FLQuant(0,dimnames=dms),
              availability   =FLQuant(1,dimnames=dms))

   ## Stock recruitment stuff
   params(res)<-FLPar(array(abPars(s = steepness, v = vbiomass, spr0 = spr0(ssb = ssb(res), rec = rec(res), fbar = fbar(res)), model=srModel),dim=c(2,1),dimnames=list(params=c("a","b"),iter=1)))

   if (is(srModel,"character")) srModel<-do.call(srModel, list())$model
   if (!is(srModel,"formula")) stop("srModel must be either formula or character")

   model(res) <-srModel
   res        <-brp(res)

   ## Equilibrium conditions at FMSY
   fbar(res)[]   <-fmsy*refpts(res)["msy","harvest",1]

   return(res)
   }

