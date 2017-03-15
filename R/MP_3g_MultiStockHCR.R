#-------------------------------------------------------------------------------
#  THIS HCR ONLY WORKS WITH SINGLE ITERATIONS!!
#
#  MultiStock HCR based on the ICES MSY framework HCR.
#
# Biomass Based HCR.
#  Reference Points: Btrigger, Blim and Fmsy.
#   No formal porposal for any of this, usually :
#           - Btrigger = Bpa
#           - Blim  = ??? Blim?YYY?
#           - Fmsy = F0.1, Fmax....
#
#  - TAC advice depending on F in relation to BRP of ALL the stocks:
#     1. First for each stock we calculate the single stock  Ftarget depending on 
#       its status in relation to the BRPs.
#           - Ftarget = Fmsy             B >= Btrigger.
#           - Ftarget = Fmsy*B/Btrigger  B <  Btrigger.
#           - Ftarget = 0.               B <  Blim (our proposal)
#     2.Calculate the ratio between Ftarget and Fsq and calculate the maximum:
#           Fadv0[st] = lambda0*Fsq : lambda0 = max(Ftarget/Fsq)
#       (=> there is only one stock for which Fadv0 = Ftarget and for the rest Fadv0 > Ftarget)

#     3.Calculate the ratio between Fupp and Fsq and calculate the minimum:
#            x_st = Fupp[st]/Fadv0[st] for all the stocks st.
#              -  If x_st >= 1 for all the stocks: lambda1 = 1 
#              -  If exist st : x_st < 1 => lambda1 = min(Fupp[st]/Fadv0[st])
#                                           and Fadv1[st] = lambda1*Fadv0 :
#     (=> there is only one stock for which Fadv0 = Ftarget and for the rest Fadv0 > Ftarget)
#
#     4.Fadv[st] = lambda1*lambda0*Fsq => TAC[st] = C[Fadv[st]]
#
#     In relation to IcesMSY HCR this new HCR has two additional arguments:
#       * advice.ctrl[['stocksInHCR']] : A vector with the name of the stocks that are taken into account in the calculation of advice.
#       * advice.ctrl[['stknm']][['ref.pts']]: A new row in the matrix with Fupp value.
#
#  ** NOTE THAT the first part of the functions does not work with iterations it should be adapted 
#        to be compatible with  multiple iterations.
#
#
#
# 08/06/2016 - Created:  Dorleta Garcia.
# 18/06/2016 - Modified: Dorleta Garcia.
# 20/12/2016 - Modified: Dorleta Garcia. Incorporation of Data Limited Stocks into the HCR.
#---------------------------------------------------------------------------------------------------
#' @rdname annualTAC
MultiStockHCR <- function(stocks, indices, advice, advice.ctrl, year, stknm,...){

   # project the stock 3 years, (current year, TAC year, TAC year + 1 for ssb or biomass constraints).
    nyears      <- ifelse(is.null(advice.ctrl[[stknm]][['nyears']]), 3, advice.ctrl[[stknm]][['nyears']])
    wts.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['wts.nyears']]), 3, advice.ctrl[[stknm]][['wts.nyears']])
    fbar.nyears <- ifelse(is.null(advice.ctrl[[stknm]][['fbar.nyears']]), 3, advice.ctrl[[stknm]][['fbar.nyears']])
    f.rescale   <- ifelse(is.null(advice.ctrl[[stknm]][['f.rescale']]), TRUE, advice.ctrl[[stknm]][['f.rescale']])
   # disc.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['disc.nyears']]), wts.nyears, advice.ctrl[[stknm]][['disc.nyears']])

    stk <- stocks[[stknm]]
    stk@harvest[stk@harvest < 0] <- 0.00001
    
    ageStruct <- ifelse(dim(stk@m)[1] > 1, TRUE, FALSE)

    stk <- FLAssess::stf(stk, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = f.rescale) #, disc.nyrs = disc.nyears)

   # if(dim(stk@m)[1] == 1)    stk@harvest[] <- stk@catch.n[]/stk@stock.n[] 
    
    ref.pts <- advice.ctrl[[stknm]]$ref.pts # matrix[6,it]  rows = Bmsy, MSY, alpha_0, alpha_1, alpha_2, beta
    Cadv <- ifelse(advice.ctrl[[stknm]][['AdvCatch']][year+1] == TRUE, 'catch', 'landings')
   
    iter     <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)

    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

   stocksInHCR    <- advice.ctrl[['stocksInHCR']]
   stocksCat      <- advice.ctrl[['stocksCategory']]
   
    # Build fwd.ctrl.
    #-----------------
   
   # Find where is the SSB (Age structured) OR Biomass (Aggregated) in relation to reference points
   # Argument!
    # Last SSB (Age structured) OR Biomass (Aggregated) estimate
    
    Ftg <- Fupp  <- Fsq <- numeric(length(stocksInHCR))
   names(Ftg) <-names(Fupp) <- names(Fsq) <- stocksInHCR
   
   
    for(st in stocksInHCR){
      
      ref.pts_st <- advice.ctrl[[st]][['ref.pts']]
      
      if(stocksCat[st] == 1){     
        if(ageStruct)
            b.datyr <- ssb(stk)[,year-1,drop = TRUE] # [it]
        else
            b.datyr <- (stk@stock.n*stk@stock.wt)[,year-1,drop = TRUE] # [it]

      
        b.pos <- apply(matrix(1:iter,1,iter),2, function(i) findInterval(b.datyr[i], ref.pts_st[c('Blim', 'Btrigger'),i]))  # [it]

        Ftg[st] <- ifelse(b.pos == 0, 0, ifelse(b.pos == 1, ref.pts['Fmsy',]*b.datyr/ref.pts_st[ 'Btrigger',], ref.pts['Fmsy',]))
        
        minfbar <- stocks[[st]]@range['minfbar']
        maxfbar <- stocks[[st]]@range['maxfbar']
        
        Fsq[st] <- mean(yearMeans(stocks[[st]]@harvest[,(year-3):(year-1)])[ac(minfbar:maxfbar),drop=T])
    
        Fupp[st] <- ref.pts_st['Fupp',]
      }
      
      if(stocksCat[st] == 3){     
             Brat    <- c(mean(indices[[st]][[1]]@index[,(year-2):(year-1)])/mean(indices[[st]][[1]]@index[,(year-3):(year-5)])) # [it]
             C       <- yearMeans(stocks[[st]]@catch[,(year-3):(year-1), drop=T])   # [it]
             tac     <- advice[['TAC']][st, year-1,drop=T]
             alpha   <- advice.ctrl[[stknm]][["ref.pts"]]["alpha", ]
             beta    <- advice.ctrl[[stknm]][["ref.pts"]]["beta", ]
             betaUp  <- advice.ctrl[[stknm]][["ref.pts"]]["betaUp", ]
             
             tacUpMult <- ifelse(Brat < 1-alpha, 1,      ifelse(Brat < 1+alpha, 1 + beta, 1 + betaUp))
             tacMult   <- ifelse(Brat < 1-alpha, 1-beta, ifelse(Brat < 1+alpha, 1, 1 + beta))
             
             # translate the TAC multiplier to C multiplier.
             CupMult <- tacUpMult*tac/C
             CMult   <- tacMult*tac/C
             
             # For this stocks we fill the Fsq, Fupp and Ftg in terms of catch, the linearity is assumed in
             # effort catch becasue we don't have any other information.
             Ftg[st]  <- CMult 
             Fupp[st] <- CupMult 
             Fsq[st]  <- C     
      }
    }
   
   
   
    
    lambda0 <- max(Ftg/Fsq)
   
    Fadv0 <- Fadv <- Fsq*lambda0
   
    if(any(Fadv0 > Fupp)){
      lambda1 <- min(Ftg/Fadv0) # The F multiplier.
      Fadv <- Fadv0*lambda1
    }
     
    Fadv_st <- Fadv[stknm]
   
    int.yr <- advice.ctrl[[stknm]]$intermediate.year

    for(i in 1:iter){
      
      if(stocksCat[st] == 1){
      
        if(is.na(Fadv_st) | Fadv_st == 0){
            advice[['TAC']][stknm,year+1,,,,i] <- 0
            next
        }

        int.yr <- ifelse(is.null(int.yr), 'Fsq', int.yr)
        
        if(int.yr == 'Fsq')
            fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1),  val = c(1, Fadv_st), quantity = c( 'f', 'f'), rel.year = c(-1,NA)))
        else
            fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1),  val = c(advice$TAC[stknm,year, drop=TRUE][i], Fadv_st), quantity = c( 'catch', 'f')))

        # Refresh the years in fwd!!
        fwd.ctrl@target$year     <- fwd.ctrl@target$year + assyrnumb
        fwd.ctrl@target$rel.year <- fwd.ctrl@target$rel.year + assyrnumb
    
        stki <- iter(stk, i)

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
                
                if(sr.model != 'geomean') sr.pars <- try(params(sr(FLSR(rec = rec, ssb = ssb, model = sr.model))), silent = TRUE) 
                
                if(class(sr.pars) == 'try-error' | sr.model == 'geomean'){
                    sr.model <- 'geomean'
                    sr.pars <- c(prod(c(rec))^(1/length(c(rec))))
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
            
            stki <- fwdBD(stki, fwd.ctrl, growth.years)
        }
     
        yy <- ifelse(slot(stki, Cadv)[,year+1] == 0, 1e-6, slot(stki, Cadv)[,year+1])
     
        advice[['TAC']][stknm,year+1,,,,i] <- yy # The TAC is given in terms of CATCH.

#        cat('---------------- HCR------------------------\n')
#        cat(c(fbar(stki)[,(year-1):year]), '\n')
#        cat('-------------------------------------------\n')
   #     save(stki, file = 'stki.RData')
      }
      if(stocksCat[st] == 3){
        advice[['TAC']][stknm,year+1,,,,i] <- Fsq[st]*Fadv_st
      }


      }
    return(advice)
}
