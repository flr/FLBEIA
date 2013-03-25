#-------------------------------------------------------------------------------
#                           MAIN FUNCTION
#
# Dorleta Garcia
# Created: 20/12/2010 21:07:45
# Changes: 
#   * 2012-06-15 12:14:20  Sonia Sánchez - for allowing assessment in different seasons and multiannual advice
#   * 12/03/2013 10:16:55  Dorleta Garcia -  Default value (NULL) for optional objects.
#-------------------------------------------------------------------------------

BEIA <- function(biols, SRs = NULL, BDs = NULL, fleets, covars = NULL, indices = NULL, advice = NULL, main.ctrl, biols.ctrl, 
        fleets.ctrl, covars.ctrl, obs.ctrl, assess.ctrl,  advice.ctrl){

    # Extract the common dimensions [year, season, it] from the 1st Biol.
    ny <- dim(biols[[1]]@n)[2]
    ns <- dim(biols[[1]]@n)[4]
    it <- dim(biols[[1]]@n)[6]
    minyear <- ac(dims(biols[[1]])$minyear)
    maxyear <- ac(dims(biols[[1]])$maxyear)
    seasons <- 1:ns
    
    # Stock names
    stnms <- names(biols)
   
    # Check that all FLQuants have the rigth [ny,ns,it] dimensions. 
    chckdim0 <- checkDims(biols,  minyear, maxyear, ns, it)
    chckdim1 <- checkDims(fleets, minyear, maxyear, ns, it)
    if(!is.null(covars)) chckdim2 <- checkDims(covars, minyear, maxyear, ns, it)
       
    # Extract years, check and convert into positions.
    sim.years <- as.numeric(main.ctrl$sim.years)
    if(!(sim.years[1] %in% as.numeric(minyear):as.numeric(maxyear))) stop('First simulation year is outside year range in the objects')
    if(!(sim.years[length(sim.years)] %in% as.numeric(minyear):as.numeric(maxyear))) stop('Last simulation year is outside year range in the objects')
    # convert sim.years in positon is the FLR objects.
    sim.years <- which(sim.years[1] == as.numeric(minyear):as.numeric(maxyear)):which(sim.years[2] == as.numeric(minyear):as.numeric(maxyear))

    
    stocks         <- vector('list', length(stnms)) 
    names(stocks) <- stnms
     
    for(yr in sim.years){
      for(ss in seasons){
            cat('############################################################\n')
            cat('-                   Year: ', yr, ', Season: ',ss, '\n')
            cat('############################################################\n')
            #~~~~~~~~~~~~~~~~ OPERATING MODELS (seasonal) ~~~~~~~~~~~~~~~~~~~~~#

        cat('************ OPERATING MODEL***************************\n')
        
        cat('------------ BIOLOGICAL OM ------------\n')
            # - Biologic OM.
            res   <- biols.om (biols = biols, fleets = fleets, SRs = SRs, BDs = BDs, covars = covars, biols.ctrl = biols.ctrl, year = yr, season = ss)
            biols <- res$biols
            SRs   <- res$SRs
            BDs   <- res$BDs

        cat('------------ FLEETS OM ------------\n')
            # - Fleets OM.
            res        <- fleets.om(fleets = fleets, biols = biols, covars = covars, advice = advice, fleets.ctrl = fleets.ctrl, advice.ctrl = advice.ctrl, year = yr, season = ss)
            fleets     <- res$fleets
            fleets.ctrl <- res$fleets.ctrl
            covars     <- res$covars

        cat('------------ COVARS OM ------------\n')
            # - Covariables OM. (the covariables can affect the covariables themselfs but also the biols and fleets, biols.ctrl and fleets.ctrl)
            res    <- covars.om(fleets = fleets, biols = biols, SRs = SRs, covars = covars, advice = advice, covars.ctrl = covars.ctrl, year = yr, season = ss)
            covars <- res$covars
            biols  <- res$biols
            fleets <- res$fleets
            SRs    <- res$SRs
        
        
        # In last year of the simulation, if last season, there is no assessment => go to the end.
        if(yr == sim.years[length(sim.years)] & ss == ns) next    
            
        for (st in stnms) {
          
          ass.yr <- advice.ctrl[[st]][['ass.year']] # assessment years
          if (is.null(ass.yr)) { # no value, then assessment yearly
              ass.yr <- sim.years
          } else if (ass.yr=='all' | is.na(ass.yr)) {
              ass.yr <- sim.years
          } else { # convert assessment years into positions
              ass.yr <- as.numeric(ass.yr)
              if(sum(!(ass.yr %in% as.numeric(minyear):as.numeric(maxyear)))>0) # check
                stop("Assessment years for: '", st, "' outside year range in the objects")
              # convert ass.yr in positon of the FLR objects
              ass.yr <- which(ass.yr == as.numeric(minyear):as.numeric(maxyear))
              for (i in 1:length(ass.yr)) ass.yr[i] <- which(ass.yr[i] == as.numeric(minyear):as.numeric(maxyear))
          }
          ass.ss <- advice.ctrl[[st]][['ass.season']]
          if (is.null(ass.ss)) { ass.ss <- ns } else if (is.na(ass.ss)) { ass.ss <- ns }
          if (!(ass.ss %in% seasons)) stop("Assessment season for: '", st, "' outside season range in the objects")
          
          if (yr %in% ass.yr & ss == ass.ss) {
            
            yr.man <- ifelse( ass.ss==ns, yr, yr+1)
      
            #~~~~~~~~~~~~~~~~ MANAGEMENT PROCEDURE.  (>=annual) ~~~~~~~~~~~~~~~#
            cat('************ MANAGEMENT PROCEDURE ****************************\n')
        
            # - Observation.
            cat('----------- OBSERVATION MODEL ------------\n')
            res          <- observation.mp(biols = biols, fleets = fleets, covars = covars, indices = indices, 
                                advice = advice, obs.ctrl = obs.ctrl, year = yr.man, season=ss, stknm=st)
            stocks[[st]] <- res$stock
            fleets.obs   <- res$fleets.obs
            indices      <- res$indices
                
            # - Assessment.
            cat('------------ ASSESSMENT MODEL ------------\n')
            datayr <- dimnames(biols[[1]]@n)[[2]][yr.man-1]
          
            stocks <- assessment.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, assess.ctrl = assess.ctrl, datayr = datayr, stknm=st)    
                
            # - Advice. 
            cat('----------------- ADVICE -----------------\n')
            advice <- advice.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, covars = covars, 
                                advice = advice, advice.ctrl = advice.ctrl, year = yr, season = ss, stknm=st)
        
          }
        }
      }
      #  browser()
    }
    
    if(!exists('stocks'))  stocks <- NULL
    
    return(list(biols = biols, fleets = fleets, covars = covars,  advice = advice, stocks = stocks, indices = indices,  fleets.ctrl = fleets.ctrl))
}
