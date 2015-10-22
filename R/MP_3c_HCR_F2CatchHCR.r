#-------------------------------------------------------------------------------
# The annual fishing mortality is given as input data and the HCR 
# just translates it into catch.
#
# The function is a copy-paste from IcesMSY HCR, but:
#
#   1. Target F is directly taken from ref.pts['Ftarget',year+1,]
#   
#
# This HCR[['Ftarget']]

# 28/03/2014 13:15:58
#-------------------------------------------------------------------------------

F2CatchHCR <- function(stocks, advice, advice.ctrl, year, stknm,...){

   # project the stock 3 years, (current year, TAC year, TAC year + 1 for ssb or biomass constraints).
    nyears      <- ifelse(is.null(advice.ctrl[[stknm]][['nyears']]), 3, advice.ctrl[[stknm]][['nyears']])
    wts.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['wts.nyears']]), 3, advice.ctrl[[stknm]][['wts.nyears']])
    fbar.nyears <- ifelse(is.null(advice.ctrl[[stknm]][['fbar.nyears']]), 3, advice.ctrl[[stknm]][['fbar.nyears']])
    f.rescale   <- ifelse(is.null(advice.ctrl[[stknm]][['f.rescale']]), TRUE, advice.ctrl[[stknm]][['f.rescale']])
   # disc.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['disc.nyears']]), wts.nyears, advice.ctrl[[stknm]][['disc.nyears']])

    stk <- stocks[[stknm]]
    stk@harvest[stk@harvest < 0] <- 0.00001
    
    ageStruct <- ifelse(dim(stk@m)[1] > 1, TRUE, FALSE)

    stk <- stf(stk, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = f.rescale) #, disc.nyrs = disc.nyears)

   # if(dim(stk@m)[1] == 1)    stk@harvest[] <- stk@catch.n[]/stk@stock.n[] 
    
    ref.pts <- advice.ctrl[[stknm]]$ref.pts # matrix[6,it]  rows = Bmsy, MSY, alpha_0, alpha_1, alpha_2, beta
    Cadv <- ifelse(advice.ctrl[[stknm]][['AdvCatch']][year+1] == TRUE, 'catch', 'landings')
   
    iter     <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)

    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

    # Build fwd.ctrl.
    #-----------------
    # Last SSB (Age structured) OR Biomass (Aggregated) estimate
    if(ageStruct)
        b.datyr <- ssb(stk)[,year-1,drop = TRUE] # [it]
    else
        b.datyr <- (stk@stock.n*stk@stock.wt)[,year-1,drop = TRUE] # [it]

    Ftg <- ref.pts['Ftarget',year+1,]
    
    print(Ftg)
    
    int.yr <- advice.ctrl[[stknm]]$intermediate.year

    for(i in 1:iter){
    
        if(is.na(Ftg[i]) | Ftg[i] == 0){
            advice[['TAC']][stknm,year+1,,,,i] <- 0
            next
        }

        int.yr <- ifelse(is.null(int.yr), 'Fsq', int.yr)
        
        if(int.yr == 'Fsq')
            fwd.ctrl <- fwdControl(data.frame(year = c(0, 1),  val = c(1, Ftg[i]), quantity = c( 'f', 'f'), rel.year = c(-1,NA)))
        else
            fwd.ctrl <- fwdControl(data.frame(year = c(0, 1),  val = c(advice$TAC[stknm,year, drop=TRUE][i], Ftg[i]), quantity = c( 'catch', 'f')))

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
                if(is.null(advice.ctrl[[stknm]]$sr$years)) sr.yrs <- which(round(quantSums(stocks[[stknm]]@stock.n))!=0)[1]:(year-1) # yr0 missing => use all data years.
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

            stki <- fwd(stki, ctrl = fwd.ctrl, sr = list(model =sr.model, params = sr1))
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
    return(advice)
}

