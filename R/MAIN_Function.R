#-------------------------------------------------------------------------------
#                           MAIN FUNCTION
#
# Dorleta Garcia
# Created: 20/12/2010 21:07:45
# Changed:17/01/2011 12:05:59
#-------------------------------------------------------------------------------

BEIA <- function(biols, SRs, BDs, fleets, covars, indices, advice, main.ctrl, biols.ctrl, 
        fleets.ctrl, covars.ctrl, obs.ctrl, assess.ctrl,  advice.ctrl){

    # Extract the coommon dimensions [year, season, it] from the 1st Biol.
    ny <- dim(biols[[1]]@n)[2]
    ns <- dim(biols[[1]]@n)[4]
    it <- dim(biols[[1]]@n)[6]
    minyear <- ac(dims(biols[[1]])$minyear)
    maxyear <- ac(dims(biols[[1]])$maxyear)
    seasons <- 1:ns
   
    # Check that all FLQuants have the rigth [ny,ns,it] dimensions. 
    chckdim0 <- checkDims(biols,  minyear, maxyear, ns, it)
    chckdim1 <- checkDims(fleets, minyear, maxyear, ns, it)
    if(!is.null(covars)) chckdim2 <- checkDims(covars, minyear, maxyear, ns, it)
       
    # Extract years, check and convert into positions.
    sim.years <- as.numeric(main.ctrl$sim.years)
    if(!(sim.years[1] %in% as.numeric(minyear):as.numeric(maxyear))) stop('First simulation year is outside year range in the objects')
    if(!(sim.years[1] %in% as.numeric(minyear):as.numeric(maxyear))) stop('Last simulation year is outside year range in the objects')
    # convert sim.years in poisiton is the FLR objects.
    sim.years <- which(sim.years[1] == as.numeric(minyear):as.numeric(maxyear)):which(sim.years[2] == as.numeric(minyear):as.numeric(maxyear))

     
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
            res        <- fleets.om(fleets = fleets, biols = biols, covars = covars, advice = advice, fleets.ctrl = fleets.ctrl, year = yr, season = ss)
            fleets     <- res$fleets
            fleets.ctrl <- res$fleets.ctrl
            covars     <- res$covars

        cat('------------ COVARS OM ------------\n')
            # - Covariables OM. (the covariables can affect the covariables themselfs but also the biols and fleets, biols.ctrl and fleets.ctrl)
            res    <- covars.om(fleets = fleets, biols = biols, covars = covars, advice = advice, covars.ctrl = covars.ctrl, year = yr, season = ss)
            covars <- res$covars
            biols  <- res$biols
            fleets <- res$fleets
        }
             
        # In last year of the simulation, there is no assessment => go to the end.
        
        if(yr == sim.years[length(sim.years)]) next
        
        #~~~~~~~~~~~~~~~~ MANAGEMENT PROCEDURE.  (>=annual) ~~~~~~~~~~~~~~~#
        cat('************ MANAGEMENT PROCEDURE ****************************\n')
    
        # - Observation.
        cat('----------- OBSERVATION MODEL ------------\n')
        res         <- observation.mp(biols = biols, fleets = fleets, covars = covars, indices = indices, 
                            advice = advice, obs.ctrl = obs.ctrl, year = yr)
        stocks      <- res$stocks
        fleets.obs  <- res$fleets.obs
        indices     <- res$indices
            
        # - Assessment.
        cat('------------ ASSESSMENT MODEL ------------\n')
        datayr <- dimnames(biols[[1]]@n)[[2]][yr-1]
      
        stocks <- assessment.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, assess.ctrl = assess.ctrl, datayr = datayr)    
            
        # - Advice. 
        cat('----------------- ADVICE -----------------\n')
        advice <- advice.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, covars = covars, 
                            advice = advice, advice.ctrl = advice.ctrl, year = yr, season = ss)
      #  browser()
    }
    
    if(!exists('stocks'))  stocks <- NULL
    
    return(list(biols = biols, fleets = fleets, covars = covars,  advice = advice, stocks = stocks, indices = indices,  fleets.ctrl = fleets.ctrl))
}










