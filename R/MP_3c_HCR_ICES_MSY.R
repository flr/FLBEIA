#-------------------------------------------------------------------------------
# HCR proposed by ICES under MSY framework
# Biomass Based HCR.
#  Reference Points: Btrigger, Blim and Fmsy.
#   No formal porposal for any of this, usually :
#           - Btrigger = Bpa
#           - Blim  = ??? Blim?YYY?
#           - Fmsy = F0.1, Fmax....
#
#  - TAC advice depending on F in relation to BRP is:
#           - TAC[Fmsy]             B >= Btrigger.
#           - TAC[Fmsy*B/Btrigger]  B <  Btrigger.
#           - 0.                    B <  Blim (our proposal)
#
# 07/09/2011 12:20:24
#-------------------------------------------------------------------------------
#' @rdname annualTAC
#' @aliases IcesHCR
IcesHCR <- function(stocks, advice, advice.ctrl, year, stknm,...){

   # project the stock 3 years, (current year, TAC year, TAC year + 1 for ssb or biomass constraints).
    nyears      <- ifelse(is.null(advice.ctrl[[stknm]][['nyears']]), 3, advice.ctrl[[stknm]][['nyears']])
    wts.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['wts.nyears']]), 3, advice.ctrl[[stknm]][['wts.nyears']])
    fbar.nyears <- ifelse(is.null(advice.ctrl[[stknm]][['fbar.nyears']]), 3, advice.ctrl[[stknm]][['fbar.nyears']])
    f.rescale   <- ifelse(is.null(advice.ctrl[[stknm]][['f.rescale']]), TRUE, advice.ctrl[[stknm]][['f.rescale']])
   # disc.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['disc.nyears']]), wts.nyears, advice.ctrl[[stknm]][['disc.nyears']])
  
    # are there constraints in TAC variability?
    TACvar <- advice.ctrl[[stknm]][['TACvar']]
    if(!is.null(TACvar)){
        if(length(TACvar) == 1) TACvar <- c(-TACvar,TACvar) # Symmetric TAC variability
        if((TACvar[1] > 0) | (TACvar[2] < 0) ) stop('TACvar must be a vector with two numbers, 
                      the first negative with the constraint when TAC decreases, and the second one 
                      positive with the constraint when TAC increases. For a symmetric constraint a single 
                      positive value can be used.')
    }

    # Fill the 0-s and NA-s with almost 0 values to avoid problems when the fishery is closed for example, or there is no catch...
    stk <- stocks[[stknm]]
    stk@harvest[stk@harvest < 1e-12 | is.na(stk@harvest)] <- 1e-12
    
    stk@catch.n[is.na(stk@catch.n)] <- 1e-6
    stk@landings.n[is.na(stk@landings.n)] <- 0
    stk@discards.n[is.na(stk@discards.n)] <- 1e-6
    
    stk@catch.n[stk@catch.n==0] <- 1e-6
    stk@landings.n[stk@landings.n==0] <- 1e-6
    stk@discards.n[stk@discards.n==0] <- 0
    
    stk@catch <- computeCatch(stk)
    stk@catch[is.na(stk@catch)] <- 0
    stk@landings <- computeLandings(stk)
    stk@discards <- computeDiscards(stk)
    
    ageStruct <- ifelse(dim(stk@m)[1] > 1, TRUE, FALSE)
    


    if(ageStruct == TRUE){
      if(any(stk@catch[,tail(dimnames(stk@catch)$year,fbar.nyears)]<1e-2)){
        stk <- stf_correctSel(stk, nyears = 3, wts.nyears = wts.nyears, fbar.nyears = fbar.nyears, f.rescale = f.rescale) #, disc.nyrs = disc.nyears)
        }else{
         stk <- stf(stk, nyears = 3, wts.nyears = wts.nyears, fbar.nyears = fbar.nyears, f.rescale = f.rescale) #, disc.nyrs = disc.nyears)
    
      }}else{
       stk <- stfBD(stk, nyears = 3, wts.nyears = wts.nyears, fbar.nyears = fbar.nyears)}
    
   # if(dim(stk@m)[1] == 1)    harvest(stk) <- stk@catch.n/stk@stock.n
    
    ref.pts <- advice.ctrl[[stknm]]$ref.pts # matrix[3,it]  rows = Blim, Btrigger, Fmsy
    Cadv <- ifelse(advice.ctrl[[stknm]][['AdvCatch']][year+1] == TRUE, 'catch', 'landings')
    
    # Target when the biomass is above MSY BTrigger, it can be Fupp or Fmsy, or other target in the reference points
    # If the argument is missing the default is Fmsy
    target <- ifelse(is.null(advice.ctrl[[stknm]]$target), 'Fmsy', advice.ctrl[[stknm]]$target) 
      
    iter     <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)

    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

    # Build fwd.ctrl.
    #-----------------
    
    int.yr <- advice.ctrl[[stknm]]$intermediate.year

    for(i in 1:iter){
    
        stki <- iter(stk, i)

        int.yr <- ifelse(is.null(int.yr), 'Fsq', int.yr)
        
        # For Ftg we first use Fmsy and then rerun fwd using the updated Ftg depending on the value of SSB
        Ftg <- ref.pts['Fmsy',i]
        
        if(int.yr == 'Fsq') {
		
		# Calculate Fsq
		if(fbar.nyears == 1 | f.rescale) {
		Fsq <- mean(fbar(stki)[,(year-1)]) 
		} else {
		Fsq <- mean(fbar(stki)[,(year-fbar.nyears):(year-1)])
		}

            fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1),  val = c(Fsq, Ftg), quantity = c( 'f', 'f'), rel.year = c(NA,NA))) 
	}   else {
            fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1),  val = c(advice$TAC[stknm,year, drop=TRUE][i], Ftg), quantity = c( 'catch', 'f')))
	}

        # Refresh the years in fwd!!
        fwd.ctrl@target$year     <- fwd.ctrl@target$year + assyrnumb
        fwd.ctrl@target$rel.year <- fwd.ctrl@target$rel.year + assyrnumb
    

        # if in <year 0> quantity = catch => set TAC in <year 0> in val
        if(fwd.ctrl@target[fwd.ctrl@target$year == assyrnumb,'quantity'] == 'catch'){
            k <- which(fwd.ctrl@target$year == assyrnumb)
            fwd.ctrl@target[k,'val']     <- advice$TAC[stknm,year,,,,i]
            fwd.ctrl@trgtArray[k, 'val',] <- advice$TAC[stknm,year,,,,i]
        }

     #   if(stknm == 'CMON') browser()

        if(dim(stki@m)[1] > 1){
            # First estimate/extract the SR model and params.
            sr.pars  <- advice.ctrl[[stknm]]$sr$params # sr parameters if specified.
            sr.model <- advice.ctrl[[stknm]]$sr$model  # sr model, mandatory.
            if(is.null(sr.pars)){                   # if params missing => estimate the parameters using the specified years.
                if(is.null(advice.ctrl[[stknm]]$sr$years)) sr.yrs <- which(round(quantSums(stocks[[stknm]]@stock.n))!=0)[1]:(year-1)# yr0 missing => use all data years, except the assessment year for which rec is unknown
                else{
                    y.rm <- as.numeric(advice.ctrl[[stknm]]$sr$years['y.rm'])
                    nyrs <- as.numeric(advice.ctrl[[stknm]]$sr$years['num.years'])
                    sr.yrs <- yrsnames[(year-y.rm-nyrs + 1):(year-y.rm)]
                }
                rec <- stki@stock.n[1,sr.yrs]
                ssb <- ssb(stki)[,sr.yrs]

                # if rec.age != 0 adjust rec and ssb.
                rec.age <- as.numeric(dimnames(rec)[[1]])
                if(rec.age != 0){
                    rec <- rec[, -(1:rec.age),]
                    ssb <- ssb[, 1:(dim(ssb)[2] - rec.age),]
                }
                
                if(sr.model != 'geomean') sr.pars <- try(params(fmle(FLSR(rec = rec, ssb = ssb, model = sr.model))), silent = TRUE) 
                
                if(class(sr.pars) == 'try-error' | sr.model == 'geomean'){
                    sr.model <- 'geomean'
                    sr.pars <- 10^6*c(prod(c(rec/10^6))^(1/length(c(rec))))
                    sr.pars <- FLPar(a = ifelse(is.na(sr.pars), 0, sr.pars))
                }
                 
                sr1 <- sr.pars
            }
            else{ # sr.pars not null
                if(i == 1){
                   sr1 <- iter(sr.pars,i)
                }
                sr1[] <-  iter(sr.pars,i)[]

            }

            stk.aux <- FLash::fwd(stki, ctrl = fwd.ctrl, sr = list(model =sr.model, params = sr1))
            
            
            # SSB in the advice year.
            b.datyr <- ssb(stk.aux)[,year+1,drop = TRUE] # [1]
            
            # Find where the SSB (Age structured) OR Biomass (Aggregated) in relation to reference points.
            b.pos <- findInterval(b.datyr, ref.pts[c('Blim', 'Btrigger'),i])  # [1]
            Ftg <- ifelse(b.pos == 0, 0, ifelse(b.pos == 1, ref.pts['Fmsy',i]*b.datyr/ref.pts[ 'Btrigger',i], ref.pts[target,i]))
            
            print(Ftg)
            
            
            if(is.na(Ftg) | Ftg == 0){
              advice[['TAC']][stknm,year+1,,,,i] <- 0
              next
            }
            
            # Update the control and rerun the projection
            fwd.ctrl@target[2,'val'] <- fwd.ctrl@trgtArray[2,'val',] <- Ftg
            stki <- FLash::fwd(stki, ctrl = fwd.ctrl, sr = list(model =sr.model, params = sr1))
            
        }
        else{

            # Extract the years to calculate the mean historical growth of the stock
            if(is.null(advice.ctrl[[stknm]]$growth.years))   growth.years <- max(1,(year - 11)):(year-1)
            else{
                y.rm <- as.numeric(advice.ctrl[[stknm]]$growth.years['y.rm'])
                nyrs  <- as.numeric(advice.ctrl[[stknm]]$growth.years['num.years'])
                growth.years <- yrsnames[(year-y.rm-nyrs + 1):(year-y.rm)]
            }
            
            stk.aux <- fwdBD(stki, fwd.ctrl, growth.years)
            
            # SSB in the advice year.
            b.datyr <- (stk.aux@stock.n*stk.aux@stock.wt)[,year+1,drop = TRUE] # [1]
            
            # Find where the SSB (Age structured) OR Biomass (Aggregated) in relation to reference points.
            b.pos <-  findInterval(b.datyr, ref.pts[c('Blim', 'Btrigger'),i])  # [1]
            
            Ftg <- ifelse(b.pos == 0, 0, ifelse(b.pos == 1, ref.pts['Fmsy',i]*b.datyr/ref.pts[ 'Btrigger',i], ref.pts[target,i]))
            
            print(Ftg)
            
            # Update the control and rerun the projection
            fwd.ctrl@target[2,'val'] <- fwd.ctrl@trgtArray[2,'val',] <- Ftg
            stki <- fwdBD(stki, fwd.ctrl, growth.years)
            
        }
     
        yy  <- ifelse(slot(stki, Cadv)[,year+1] == 0, 1e-6, slot(stki, Cadv)[,year+1])
        yy0 <- advice[['TAC']][stknm,year,,,,i]
        
        if(!is.null(TACvar) & b.pos == 2 & yy0 > 1e-6){ # If there is constraint, TAC[y-1] > 0 and SSB > MSY Btrigger
          int <- findInterval(yy/yy0-1, TACvar) 
          yy <- ifelse(int == 1, yy, 
                    ifelse(int == 0, yy0*(1+TACvar[1]), 
                                     yy0*(1+TACvar[2])))
        }

        advice[['TAC']][stknm,year+1,,,,i] <- yy # The TAC is given in terms of CATCH.
        
#        cat('---------------- HCR------------------------\n')
#        cat(c(fbar(stki)[,(year-1):year]), '\n')
#        cat('-------------------------------------------\n')
   #     save(stki, file = 'stki.RData')

    }
    return(advice)
}

