#-------------------------------------------------------------------------------
#
#' Run the FLBEIA bio-economic simulation model
#' 
#' FLBEIA is a simulation model that describes a fishery system under a Management Strategy Estrategy framework.The objective of the model is to facilitate the Bio-Economic evaluation of Management strategies. The model is multistock, multifleet and seasonal. The simulation is divided in 2 main blocks, the Operating Model (OM) and the Management Procedure  (MP). In turn, the OM is divided in 3 components, the biological, the fleets and the covariables component. The MP is also divided in 3 components, the observation, the assessment and the advice.
#'
#' @param biols An FLBiols object.
#' @param fleets An FLFleetsExt object. An extended version of the FLFleet object defined in FLCore. 
#' @param SRs A list of FLSRSim objects. One per age structured stock in biols object.
#' @param BDs A list of FLSRSim objects. One per biomass dynamic stock in biols object.
#' @param covars A list of FLQuants used to store any kind of variables that are used within the simulation and are not stored in the standard objects. 
#' @param indices A list of FLIndices. Each element must correspond with one of the stocks in biols object. 
#' @param advice A list with two FLQuant elements, TAC and quota.share. TAC is an FLQuant  with quant dimension equal to the number of stocks in biols object, the names used in in the quant dimension must be equal to those used in biols. quota.share is a list with one element per stock in biols object indicating the quota share per stock and fleet. The quant dimension of the elements must be equal to the number of fleets and the names used must be equal to those in fleets objects.  
#' @param main.ctrl A list with the settings to control the main function (the year range,...).
#' @param biols.ctrl A list with the settings to control the biological operating model for each stock (the population dynamic model used,  additional parameters,...)
#' @param fleets.ctrl A list with the settings to control the fleets operating model for each fleet (the fleets' short and long term dynamic models used, price model, additional parameters,...)
#' @param covars.ctrl A list with the settings to control the covars operating model for each fleet (a dynamic model for each covariable,  additional parameters, ...)
#' @param obs.ctrl A list with the settings to control the observation model for each stock (the observation model for the stock, for stock dependent indices, additional parameters, ...)
#' @param assess.ctrl A list with the settings to control the specify the assessment model for each stock (the assessment model for the stock and the control parameters used to run the model)
#' @param advice.ctrl A list with the settings to control the advice model for each stock (the HCR for each stock, the reference points used in the HCR, additional parameters, ...)
#
#' @return A list with 8 elements biols, fleets, covars,  advice, stocks, indices,  fleets.ctrl,
#'          pkgs.versions. All the elements except stocks and pkgs.versions correspond with the 
#'          the updated versions of the objects used in the call to FLBEIA. stocks is a list of FLStocks object
#'          containing the perceived stocks used in the management procedure to produce the management advice. 
#'          pkgs.versions is a matrix indicating the packages and package version used along the simulation. 
#
#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(FLAssess)          # required to use the IcesHCR. Not available for win64
#' library(FLash)             # required to use the IcesHCR. Not available for win64
#' library(ggplot2)  
#' 
#' 
#' #---------------------------------------------------------------- 
#' # Example with 1 stock, 1 Fleets, 1 seasons and 1 iteration: one
#' #----------------------------------------------------------------
#' 
#' # Load the data to run FLBEIA in a one stock one fleet example using the HCR used by ICES in the MSY framework. 
#'  data(one) 
#' 
#' # The names and the class of the objects needed to run FLBEIA.
#' # sapply(ls(), function(x) class(get(x)))
#' 
#' # In this scenario a single, age-structured, stock is exploited by a single fleet with a unique metier. 
#' # The fleet catches yearly exactly the adviced TAC and there is no exit-entry of vessels in the fishery.  
#' # The stock abundance and exploitation level is observed without error in the observation model.
#' # There is no assessment model and the TAC advice is used through the HCR used by ICES in the MSY framework.  
#'  
#' 
#'   
#'     
#' s0 <- FLBEIA(biols = oneBio,       # FLBiols object with one FLBiol element for stk1.
#'                SRs = oneSR,        # A list with one FLSRSim object for stk1.
#'                BDs = NULL,         # No Biomass Dynamic populations in this case.
#'             fleets = oneFl,        # FLFleets object with on fleet.
#'             covars = oneCv,         # covars not used
#'            indices = NULL,         # indices not used 
#'             advice = oneAdv,       # A list with two elements 'TAC' and 'quota.share'
#'          main.ctrl = oneMainC,     # A list with one element to define the start and end of the simulation.
#'         biols.ctrl = oneBioC,      # A list with one element to select the model to simulate the stock dynamics.
#'        fleets.ctrl = oneFlC,       # A list with several elements to select fleet dynamic models and store additional parameters.
#'        covars.ctrl = oneCvC,         # covars control not used 
#'           obs.ctrl = oneObsC,      # A list with one element to define how the stock observed ("PerfectObs").
#'        assess.ctrl = oneAssC,      # A list with one element to define how the stock assessment model used ("NoAssessment").
#'        advice.ctrl = oneAdvC)       # A list with one element to define how the TAC advice is obtained ("IcesHCR").
#' 
#' # Names of the object returned by FLBEIA
#' names(s0)
#' 
#' # The default plot for FLBiol defined in FLCore
#' plot(s0$biols[[1]])
#' 
#' # Create summary data frames (biological, economic, and catch)
#' proj.yr     <- 2013 
#' s0_sum      <- bioSum(s0)
#' s0_flt      <- fltSum(s0)
#' s0_fltStk   <- fltStkSum(s0)
#'
#'
#' # Create several plots and save them in the working directory using 'pdf' format and 
#' # 's0' suffix in the name.
#' 
#' 
#' plotFLBiols(s0$biols, pdfnm='s0')
#' plotFLFleets(s0$fleets, pdfnm='s0')
#' plotEco(s0, pdfnm='s0')
#' plotfltStkSum(s0, pdfnm='s0')
#' 
#' #------------------------------------------------------------ 
#' # Example with several iterations: oneIters
#' #------------------------------------------------------------
#'  
#'  # Load the same data set as before but with 3 iterations.
#'  # Run FLBEIA and plot the results
#'  
#' data(oneIt)
#'  
#' s1 <- FLBEIA(biols = oneItBio,       # FLBiols object with one FLBiol element for stk1.
#'                SRs = oneItSR,        # A list with one FLSRSim object for stk1.
#'                BDs = NULL,         # No Biomass Dynamic populations in this case.
#'             fleets = oneItFl,        # FLFleets object with on fleet.
#'             covars = oneItCv,         # covars not used
#'            indices = NULL,         # indices not used 
#'             advice = oneItAdv,       # A list with two elements 'TAC' and 'quota.share'
#'          main.ctrl = oneItMainC,     # A list with one element to define the start and end of the simulation.
#'         biols.ctrl = oneItBioC,      # A list with one element to select the model to simulate the stock dynamics.
#'        fleets.ctrl = oneItFlC,       # A list with several elements to select fleet dynamic models and store additional parameters.
#'        covars.ctrl = oneItCvC,         # covars control not used 
#'           obs.ctrl = oneItObsC,      # A list with one element to define how the stock observed ("PerfectObs").
#'        assess.ctrl = oneItAssC,      # A list with one element to define how the stock assessment model used ("NoAssessment").
#'        advice.ctrl = oneItAdvC)       # A list with one element to define how the TAC advice is obtained ("IcesHCR").
#' 
#' # Names of the object returned by FLBEIA
#' names(s1)
#' 
#' # The default plot for FLBiol defined in FLCore
#' plot(s1$biols[[1]])
#' 
#' # Create summary data frames (biological, economic, and catch)
#' proj.yr     <- 2013 
#' s1_bio     <- bioSum(s1)
#' s1_flt     <- fltSum(s1)
#' s1_fltStk  <- fltStkSum(s1)
#' 
#' s1_bioQ    <- bioSumQ(s1_bio)
#' s1_fltQ    <- fltSumQ(s1_flt)
#' s1_fltStkQ <- fltStkSumQ(s1_fltStk)
#' 
#' s1b_bio     <- bioSum(s1, long = FALSE)
#' s1b_flt     <- fltSum(s1, long = FALSE)
#' s1b_fltStk  <- fltStkSum(s1, long = FALSE)
#' 
#' s1b_fltQ    <- bioSumQ(s1b_bio)
#' s1b_fltQ    <- fltSumQ(s1b_flt)
#' s1b_fltStkQ <- fltStkSumQ(s1b_fltStk)
#' 
#' # Create several plots and save them in the working directory using 'pdf' format and 
#' # 's1' suffix in the name.
#' 
#' #' plotFLBiols(s1$biols, pdfnm='s1')
#' plotFLFleets(s1$fleets, pdfnm='s1')
#' plotEco(s1, pdfnm='s1')
#' plotfltStkSum(s1, pdfnm='s1') 
#'
#'  
#' #------------------------------------------------------------------ 
#' # Example with 2 stock, 2 Fleets, 4 seasons and 1 iteration: multi
#' #------------------------------------------------------------------
#'  
#'  # Load the multi data set. This dataset has 2 stocks, one stk1 is 
#'  # age structured and the second one stk2 is aggregated in biomass.
#'  
#' data(multi)
#'  
#'  # Run FLBEIA.
#'  
#' s2 <- FLBEIA(biols = multiBio,       # FLBiols object with 2 FLBiol element for stk1.
#'                SRs = multiSR,        # A list with 1 FLSRSim object for stk1.
#'                BDs = multiBD,        # A list with 1 FLBDSim object for stk2.
#'             fleets = multiFl,        # FLFleets object with on fleet.
#'             covars = multiCv,         # covars not used
#'            indices = NULL,         # indices not used 
#'             advice = multiAdv,       # A list with two elements 'TAC' and 'quota.share'
#'          main.ctrl = multiMainC,     # A list with one element to define the start and end of the simulation.
#'         biols.ctrl = multiBioC,      # A list with one element to select the model to simulate the stock dynamics.
#'        fleets.ctrl = multiFlC,       # A list with several elements to select fleet dynamic models and store additional parameters.
#'        covars.ctrl = multiCvC,         # covars control not used 
#'           obs.ctrl = multiObsC,      # A list with one element to define how the stock observed ("PerfectObs").
#'        assess.ctrl = multiAssC,      # A list with one element to define how the stock assessment model used ("NoAssessment").
#'        advice.ctrl = multiAdvC)       # A list with one element to define how the TAC advice is obtained ("IcesHCR").
#' 
#' # Names of the object returned by FLBEIA
#' names(s2)
#' 
#' # The default plot for FLBiol defined in FLCore
#' plot(s2$biols[[1]])
#' 
#' # Create summary data frames (biological, economic, and catch)
#' 
#' s2_sum      <- bioSum(s2)
#' s2_flt      <- fltSum(s2)
#' 
#' s2b_flt     <- fltSum(s2, byyear = FALSE)
#' 
#' s2_fltStk   <- fltStkSum(s2)
#'
#' # Create several plots and save them in the working directory using 'pdf' format and 
#' # 's2' suffix in the name.
#' 
#' plotFLBiols(s2$biols, pdfnm='s2')
#' plotFLFleets(s2$fleets, pdfnm='s2')
#' plotEco(s2, pdfnm='s2')
#' plotfltStkSum(s2, pdfnm='s2')

#' }
#' 
#' 
# MAIN FUNCTION
#
# Dorleta Garcia
# Created: 20/12/2010 21:07:45
# Changes: 
#   * 2012-06-15 12:14:20  Sonia Sanchez - for allowing assessment in different seasons and multiannual advice
#   * 12/03/2013 10:16:55  Dorleta Garcia -  Default value (NULL) for optional objects.
#   * 20/10/2016           Itsaso Carmona - Check if some arguments are missing
#-------------------------------------------------------------------------------

# @export
FLBEIA <- function(biols, SRs = NULL, BDs = NULL, fleets, covars = NULL, indices = NULL, advice = NULL, main.ctrl, biols.ctrl, 
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
    
    # If SimultaneousMngt argument missing in main.ctrl => set it to FALSE, the original FLBEIA configuration.
    main.ctrl$SimultaneousMngt <- ifelse(is.null(main.ctrl$SimultaneousMngt), FALSE, TRUE)
   
    # Check that all FLQuants have the rigth [ny,ns,it] dimensions. 
    chckdim0 <- checkDims(biols,  minyear, maxyear, ns, it)
    chckdim1 <- checkDims(fleets, minyear, maxyear, ns, it)
    if(!is.null(covars)) chckdim2 <- checkDims(covars, minyear, maxyear, ns, it)
    # Check when the model to describe BD is Pellatom, that alpha has the right values.

    if(!is.null(BDs)){
      BDnms<- names(BDs)
      for(stk.bd in BDnms){
        if(BDs[[stk.bd]]@model=="PellaTom"){
          p <- BDs[[stk.bd]]@params["p",,,]
          r <- BDs[[stk.bd]]@params["r",,,]
          K <- BDs[[stk.bd]]@params["K",,,]
          if(any(BDs[[stk.bd]]@alpha<1) || any(as.vector(BDs[[stk.bd]]@alpha) > as.vector(((p/r+1)^(1/p))))){
            stop("alpha<1 or alpha > min((p/r+1)^(1/p))")
          }}}}
    
    # Extract years, check and convert into positions.
    sim.years <- as.numeric(main.ctrl$sim.years)
    if(!(sim.years[1] %in% as.numeric(minyear):as.numeric(maxyear))) stop('First simulation year is outside year range in the objects')
    if(!(sim.years[length(sim.years)] %in% as.numeric(minyear):as.numeric(maxyear))) stop('Last simulation year is outside year range in the objects')
    # convert sim.years in positon is the FLR objects.
    sim.years <- which(sim.years[1] == as.numeric(minyear):as.numeric(maxyear)):which(sim.years[2] == as.numeric(minyear):as.numeric(maxyear))
    
    # Check if the argument LandObl is missing for any fleet in fleets.ctrl. 
    # No Landing Obligation if the argument is missing.
    for (flnm in names(fleets)){ 
      if(is.null(fleets.ctrl[[flnm]]$LandObl)) fleets.ctrl[[flnm]]$LandObl<-FALSE
    }
    # Check if the argument AdvCatch is missing for any stock in advice.ctrl. 
    # AdvCatch = FALSE (TAC in terms of landings for each year) in case AdvCatch is missing.
    for (st in names(biols)){ 
      if(is.null(advice.ctrl[[st]]$AdvCatch)){
         advice.ctrl[[st]]$AdvCatch <- rep(FALSE,ny)
         names(advice.ctrl[[st]]$AdvCatch) <- c(as.numeric(minyear):as.numeric(maxyear))
      }
    }
    
    stocks         <- vector('list', length(stnms)) 
    names(stocks) <- stnms
    
    if(main.ctrl$SimultaneousMngt == FALSE) {
      
      # Define assessment conditions:
      ass.yr <- ass.ss <- vector('list', length(stnms))
      names(ass.yr) <- names(ass.ss) <- stnms
      for (st in stnms) {
        
        # Assessment years
        ass.yr[[st]] <- advice.ctrl[[st]][['ass.year']] # assessment years
        if (is.null(ass.yr[[st]])) { # no value, then assessment yearly
          ass.yr[[st]] <- sim.years
        } else if (ass.yr[[st]]=='all' | is.na(ass.yr[[st]])) {
          ass.yr[[st]] <- sim.years
        } else { # convert assessment years into positions
          ass.yr[[st]] <- as.numeric(ass.yr[[st]])
          if(sum(!(ass.yr[[st]] %in% as.numeric(minyear):as.numeric(maxyear)))>0) # check
            stop("Assessment years for: '", st, "' outside year range in the objects")
          # convert ass.yr[[st]] in positon of the FLR objects
          for (i in 1:length(ass.yr[[st]])) ass.yr[[st]][i] <- which(ass.yr[[st]][i] == as.numeric(minyear):as.numeric(maxyear))
        }
        
        # Assessment seasons
        ass.ss[[st]] <- advice.ctrl[[st]][['ass.season']]
        if (is.null(ass.ss[[st]])) { ass.ss[[st]] <- ns } else if (is.na(ass.ss[[st]])) { ass.ss[[st]] <- ns }
        if (!(ass.ss[[st]] %in% seasons)) stop("Assessment season for: '", st, "' outside season range in the objects")
        
        # Assessment year estimates necessary?
        acy <- advice.ctrl[[st]]$ass.curryr # TRUE if estimates also for assessment year are needed
        if (is.null(advice.ctrl[[st]]$ass.curryr)) { acy <- F } else if (is.na(advice.ctrl[[st]]$ass.curryr)) { acy <- F }
        obs.ctrl[[st]]$obs.curryr <- assess.ctrl[[st]]$ass.curryr <- acy
        
      }
      
    } else { # if main.ctrl$SimultaneousMngt == TRUE:
      
      # Assessment years NOT to be defined (are assumed to be all years)
      if (!is.null(advice.ctrl[[st]][['ass.year']]))
        stop("Assessment years for: '", st, "' should not be defined if  main.ctrl$SimultaneousMngt == TRUE.
              See advice.ctrl[['", st, "']][['ass.year']].")
     
      # Assessment seasons NOT to be defined (are assumed to be only once, in the last season)
      if (!is.null(advice.ctrl[[st]][['ass.season']]))
        stop("Assessment seasons for: '", st, "' should not be defined if main.ctrl$SimultaneousMngt == TRUE. 
              See advice.ctrl[['", st, "']][['ass.season']].")

      # No possibility of doing the assessment in current year
      obs.ctrl <- lapply(obs.ctrl, function(x){
                                        x[['obs.curryr']] <- FALSE
                                        return(x)})
      assess.ctrl <- lapply(assess.ctrl, function(x){
        x[['ass.curryr']] <- FALSE
        return(x)})
      
      }
     
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
            res        <- fleets.om(fleets = fleets, biols = biols, BDs = BDs, covars = covars, advice = advice, biols.ctrl = biols.ctrl, fleets.ctrl = fleets.ctrl, advice.ctrl = advice.ctrl, year = yr, season = ss)
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
         
        if(main.ctrl$SimultaneousMngt == FALSE){   
          for (st in stnms) {
          
            if (yr %in% ass.yr[[st]] & ss == ass.ss[[st]]) {
              
              yr.man <- ifelse( ass.ss[[st]]==ns, yr, yr+1)
        
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
            
              res <- assessment.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, covars=covars, 
                                      assess.ctrl = assess.ctrl, datayr = datayr, season=ss, stknm=st)  
              stocks <- res$stocks
              covars <- res$covars
  
              # - Advice. 
              cat('----------------- ADVICE -----------------\n')
              advice <- advice.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, covars = covars, 
                                  advice = advice, advice.ctrl = advice.ctrl, year = yr, season = ss, stknm=st)
      
        }}}}
        if(main.ctrl$SimultaneousMngt == TRUE & yr < sim.years[length(sim.years)]){  # Simultaneous and Yearly management. 
        
        #~~~~~~~~~~~~~~~~ MANAGEMENT PROCEDURE.  (>=annual) ~~~~~~~~~~~~~~~#
        cat('************ MANAGEMENT PROCEDURE ****************************\n')
        
        stocks <- vector('list', length(stnms))
        names(stocks) <- stnms
        #indices <- vector('list', length(stnms))
        #names(indices) <- stnms
        
        # - Observation.
        cat('----------- OBSERVATION MODEL ------------\n')
        for(st in stnms){

          res          <- observation.mp(biols = biols, fleets = fleets, covars = covars, indices = indices, 
                                       advice = advice, obs.ctrl = obs.ctrl, year = yr, season=ns, stknm=st)
          stocks[[st]] <- res$stock
          fleets.obs   <- res$fleets.obs
          indices      <- res$indices
        }
        
        # - Assessment.
        cat('------------ ASSESSMENT MODEL ------------\n')
        for(st in stnms){
        
          res <- assessment.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, covars = covars, 
                               assess.ctrl = assess.ctrl, datayr = yr-1, season=ns, stknm=st)    
          stocks[[st]] <- res$stocks[[st]]
          covars <- res$covars
          }
        
        
        # - Advice. 
        cat('----------------- ADVICE -----------------\n')
        for(st in stnms){
          advice <- advice.mp(stocks = stocks, fleets.obs = fleets.obs, indices = indices, covars = covars, 
                            advice = advice, advice.ctrl = advice.ctrl, year = yr, season = ns, stknm=st)
          }
        }

    }
    
    if(!exists('stocks'))  stocks <- NULL
    
    return(list(biols = biols, fleets = fleets, covars = covars,  advice = advice, stocks = stocks, indices = indices,  fleets.ctrl = fleets.ctrl,
                      pkgs.versions = installed.packages(fields = 'Packaged')[,c('Built', 'Version', 'Packaged')]))
}
