#-------------------------------------------------------------------------------
#                          HCRs
#   - MAPHCR a HCR specially designed to fulfill the requirement of Multi-Anual Managemen plans for Woth Western Waters. 
#
# Dorleta GarcYYYa LO HEMOS DEJADO EN LA LINEA 104!!!!!!!!!!
# Created: 20/04/2015 13:26:13
# Changed: 20/04/2015 13:26:18
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# MAPHCR(stocks, covars, advice, hcr.ctrl)   
# year = Assessment year (POSITION).
#
# THIS HCR ONLY WORKS FOR AGE STRUCTURED STOCKS!!

# advice.ctrl elements:
# ref.pts :: Bpa and Ftg (Fixed).
# N: The number of years to recover ssb
#-------------------------------------------------------------------------------

MAPHCR <- function(stocks, advice, advice.ctrl, year, stknm,...){
   
   # project the stock 3 years, (current year, TAC year, TAC year + 1 for ssb or biomass constraints). 
#    nyears      <- ifelse(is.null(advice.ctrl[[stknm]][['nyears']]), 3, advice.ctrl[[stknm]][['nyears']]) 
    wts.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['wts.nyears']]), 3, advice.ctrl[[stknm]][['wts.nyears']]) 
    fbar.nyears <- ifelse(is.null(advice.ctrl[[stknm]][['fbar.nyears']]), 3, advice.ctrl[[stknm]][['fbar.nyears']]) 
    f.rescale   <- ifelse(is.null(advice.ctrl[[stknm]][['f.rescale']]), TRUE, advice.ctrl[[stknm]][['f.rescale']]) 
   # disc.nyears  <- ifelse(is.null(advice.ctrl[[stknm]][['disc.nyears']]), wts.nyears, advice.ctrl[[stknm]][['disc.nyears']]) 
    
    Bpa <- advice.ctrl[[stknm]][['ref.pts']]['Bpa',year+1,]
    Ftg <- advice.ctrl[[stknm]][['ref.pts']]['Ftarget',year+1,]
    Cup <- advice.ctrl[[stknm]][['ref.pts']]['Cup',year+1,]
    Clo <- advice.ctrl[[stknm]][['ref.pts']]['Clo',year+1,]
    N <- advice.ctrl[[stknm]][['N']]
    Cadv <- ifelse(advice.ctrl[[stknm]][['AdvCatch']][year+1] == TRUE, 'catch', 'landings')
   
    stk <- stocks[[stknm]]
    stk@harvest[stk@harvest < 0] <- 0.00001
  
    stk <- stf(stk, nyears = 3, wts.nyears = wts.nyears, fbar.nyears = fbar.nyears, f.rescale = f.rescale) #, disc.nyrs = disc.nyears)
    
    iter   <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)
    
    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

    TACvar <- FALSE
   
    TACnow <- advice$TAC[stknm, year, drop=T]
   
    
   # The projection is done in two steps, first the stock is projected up to year 1st january TAC year, in order to calculate SSB 
   # in 1st January data year.
   # If this SSB > Bpa => Fobj = Ftarget.
   # else if Fadv tq SSB[tac year + 1] > SSB[data year] + (Bpa - SSB[data year])/K, if this leads to and F > Ftarget => Fadv = Farget

   
    for(i in 1:iter){
        
      stki <- iter(stk, i)
        
      # CALCULATE RECRUITMENT MODEL IN FOR THE PROJECTION
                                      
      # First estimate/extract the SR model and params.
      sr.pars  <- advice.ctrl[[stknm]]$sr$params # sr parameters if specified.  
      sr.model <- advice.ctrl[[stknm]]$sr$model  # sr model, mandatory.  
      if(is.null(sr.pars)){                   # if params missing => estimate the parameters using the specified years.
      if(is.null(advice.ctrl[[stknm]]$sr$years)){  # We remove the first year in case they are zero because it means the stock didn't exist in this year.
          sr.yrs <-which(round(quantSums(stocks[[stknm]]@stock.n))!=0)[1]:(year-1) # yr0 missing => use all data years.
      }
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
      } else { # sr.pars not null
        if(i == 1){
          sr1 <- iter(sr.pars,i)
        }
          sr1[] <-  iter(sr.pars,i)[]
         
      }
            
            
      # Project the population one year using Fsq.
      fsq <- mean(fbar(stki)[,ac((assyrnumb-3):(assyrnumb-1))])
      fwd.ctrl1 <- fwdControl(data.frame(year = assyrnumb,  val = c(fsq), quantity = c( 'f'),
                                               min = c(NA)
                                         , max  = c(NA)))
      stki <- fwd(stki, ctrl = fwd.ctrl1, sr = list(model =sr.model, params = sr1))
            
            
      # If SSB < Bpa Project the population ensuring that the biomass increases
      # SSB the biomass in the 1st January TAC year.
      ssbTACyr <- ssb(stki)[,ac(assyrnumb+1),]            
      
      if(ssbTACyr < Bpa){
        
          Above_Bpa_yr <- max(which(ssb(stki)[,ac(yrsnumbs[1]:(assyrnumb)),] > Bpa)) # Last year with SSB above Bpa.  
          Nyr_Below_Bpa <- length(yrsnumbs[(Above_Bpa_yr+1)]:(assyrnumb)) - 1 # Number of consecutive years with SSB < Bpa -1 (-1 current year outside the sum)
          
          K <-ifelse((N - Nyr_Below_Bpa) < 0, 1, (N - Nyr_Below_Bpa)) 
                 
          ssbobj <- ssbTACyr + (Bpa - ssbTACyr)/K
            
          fwd.ctrl2 <- fwdControl(data.frame(year = c(assyrnumb+1, assyrnumb+1, assyrnumb+1),  val = c(ssbobj,NA,NA), quantity = c( 'ssb', 'f', Cadv),
                                               min = c(NA, 0, TACnow[i]*Clo), max  = c(NA, Ftg, TACnow[i]*Cup)))
      
          stki <- fwd(stki, ctrl = fwd.ctrl2, sr = list(model =sr.model, params = sr1))     
      }
      { # Advice in Ftg and imposing -15% restriction.
          fwd.ctrl2 <- fwdControl(data.frame(year = c(assyrnumb+1, assyrnumb+1),  val = c(Ftg,NA), quantity = c( 'f', Cadv),
                                           min = c(NA, TACnow[i]*Clo), max  = c(NA, TACnow[i]*Cup)))
        
          stki <- fwd(stki, ctrl = fwd.ctrl2, sr = list(model =sr.model, params = sr1))           
      }
    
      yy <- ifelse(slot(stki, Cadv)[,year+1] == 0, 1e-6, slot(stki, Cadv)[,year+1])
                   
      advice[['TAC']][stknm,year+1,,,,i] <- yy # The TAC is given in terms of CATCH.

    } # Iter loop 

   
    return(advice)   
}




