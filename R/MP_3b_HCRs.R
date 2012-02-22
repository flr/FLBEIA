#-------------------------------------------------------------------------------
#                          HCRs
#   - annualTAC.
#
# Dorleta García
# Created: 20/12/2010 13:26:13
# Changed: 20/12/2010 13:26:18
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# annualTAC(stocks, covars, advice, hcr.ctrl)   
# year = Assessment year (POSITION).
#
# The targets and constraints will differ iteration by iteration, thus 'fwd' 
# must be applied iter by iter.
#-------------------------------------------------------------------------------

annualTAC <- function(stocks, advice, advice.ctrl, year, stknm,...){
   
   # project the stock 3 years, (current year, TAC year, TAC year + 1 for ssb or biomass constraints). 
    nyears      <- ifelse(is.null(advice.ctrl[[stknm]][['nyears']]), 3, advice.ctrl[[stknm]][['nyears']]) 
    wts.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['wts.nyears']]), 3, advice.ctrl[[stknm]][['wts.nyears']]) 
    fbar.nyears <- ifelse(is.null(advice.ctrl[[stknm]][['fbar.nyears']]), 3, advice.ctrl[[stknm]][['fbar.nyears']]) 
    f.rescale   <- ifelse(is.null(advice.ctrl[[stknm]][['f.rescale']]), TRUE, advice.ctrl[[stknm]][['f.rescale']]) 
   # disc.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['disc.nyears']]), wts.nyears, advice.ctrl[[stknm]][['disc.nyears']]) 
    
    stk <- stocks[[stknm]]
    stk@harvest[stk@harvest < 0] <- 0.00001
  
    stk <- stf(stk, nyears = 3, wts.nyears = 3, fbar.nyears = 3, f.rescale = TRUE) #, disc.nyrs = disc.nyears)
    
    fwd.ctrl <- advice.ctrl[[stknm]]$fwd.ctrl
    
    iter   <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)
    
    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

    # Refresh the years in fwd!! 
    fwd.ctrl@target$year     <- fwd.ctrl@target$year + assyrnumb
    fwd.ctrl@target$rel.year <- fwd.ctrl@target$rel.year + assyrnumb 
    
    for(i in 1:iter){
        
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
                if(is.null(advice.ctrl[[stknm]]$sr$years)) sr.yrs <- 1:year # yr0 missing => use all data years.
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
                sr.pars <- params(sr(FLSR(rec = rec, ssb = ssb, model = sr.model)))
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
    
        sl <- advice.ctrl[[stknm]][['advice']] # catch or landings?¿¿?¿
        advice[['TAC']][stknm,year+1,,,,i] <- slot(stki, sl)[,year+1]
        
#        cat('---------------- HCR------------------------\n')
#        cat(c(fbar(stki)[,(year-1):year]), '\n')
#        cat('-------------------------------------------\n')
   #     save(stki, file = 'stki.RData')
    
    }
    return(advice)   
}




