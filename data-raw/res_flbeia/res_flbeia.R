###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:  res_flbeia dataset (outputs of FLBEIA runs)
#                                                                 05/07/2018
###############################################################################

#==============================================================================
# Install Packages
#==============================================================================

  library(FLCore) 
  library(FLAssess)
  library(FLash)
  library(FLFleet)
  library(FLXSA)
  library(FLBEIA) 


#==============================================================================
#  one dataset                                                            ----
#==============================================================================
  
  load( "../../data/one.RData")
  
  s0 <- FLBEIA(biols = oneBio,    # FLBiols object with one FLBiol element for stk1.
               SRs = oneSR,       # A list with one FLSRSim object for stk1.
               BDs = NULL,        # No Biomass Dynamic populations in this case.
               fleets = oneFl,    # FLFleets object with on fleet.
               covars = oneCv,    # covars not used
               indices = NULL,    # indices not used 
               advice = oneAdv,   # A list with two elements 'TAC' and 'quota.share'
               main.ctrl = oneMainC,   # A list with one element to define the start and end of the simulation.
               biols.ctrl = oneBioC,   # A list with one element to select the model to simulate the stock dynamics.
               fleets.ctrl = oneFlC,   # A list with several elements to select fleet dynamic models and store additional parameters.
               covars.ctrl = oneCvC,   # covars control not used 
               obs.ctrl = oneObsC,     # A list with one element to define how the stock observed ("PerfectObs").
               assess.ctrl = oneAssC,  # A list with one element to define how the stock assessment model used ("NoAssessment").
               advice.ctrl = oneAdvC)  # A list with one element to define how the TAC advice is obtained ("IcesHCR").

  oneRes <- s0
  
  
#==============================================================================
#  oneIt dataset                                                          ----
#==============================================================================

  load( "../../data/oneIt.RData")
  
  s1 <- FLBEIA(biols = oneItBio,   # FLBiols object with oneIt FLBiol element for stk1.
               SRs = oneItSR,      # A list with oneIt FLSRSim object for stk1.
               BDs = NULL,         # No Biomass Dynamic populations in this case.
               fleets = oneItFl,   # FLFleets object with on fleet.
               covars = oneItCv,   # covars not used
               indices = NULL,     # indices not used 
               advice = oneItAdv,  # A list with two elements 'TAC' and 'quota.share'
               main.ctrl = oneItMainC,  # A list with oneIt element to define the start and end of the simulation.
               biols.ctrl = oneItBioC,  # A list with oneIt element to select the model to simulate the stock dynamics.
               fleets.ctrl = oneItFlC,  # A list with several elements to select fleet dynamic models and store additional parameters.
               covars.ctrl = oneItCvC,  # covars control not used 
               obs.ctrl = oneItObsC,    # A list with oneIt element to define how the stock observed ("PerfectObs").
               assess.ctrl = oneItAssC, # A list with oneIt element to define how the stock assessment model used ("NoAssessment").
               advice.ctrl = oneItAdvC) # A list with oneIt element to define how the TAC advice is obtained ("IcesHCR").
  
  oneItRes <- s1
  
  
#==============================================================================
#  multi dataset                                                            ----
#==============================================================================
  
  load( "../../data/multi.RData")
  
  s2 <- FLBEIA(biols = multiBio,   # FLBiols object with 2 FLBiol element for stk1.
               SRs = multiSR,      # A list with 1 FLSRSim object for stk1.
               BDs = multiBD,      # A list with 1 FLBDSim object for stk2.
               fleets = multiFl,   # FLFleets object with on fleet.
               covars = multiCv,   # covars not used
               indices = NULL,     # indices not used 
               advice = multiAdv,  # A list with two elements 'TAC' and 'quota.share'
               main.ctrl = multiMainC,  # A list with one element to define the start and end of the simulation.
               biols.ctrl = multiBioC,  # A list with one element to select the model to simulate the stock dynamics.
               fleets.ctrl = multiFlC,  # A list with several elements to select fleet dynamic models and store additional parameters.
               covars.ctrl = multiCvC,  # covars control not used 
               obs.ctrl = multiObsC,    # A list with one element to define how the stock observed ("PerfectObs").
               assess.ctrl = multiAssC, # A list with one element to define how the stock assessment model used ("NoAssessment").
               advice.ctrl = multiAdvC) # A list with one element to define how the TAC advice is obtained ("IcesHCR").
  
  multiRes <- s2
 
  
#==============================================================================
# Save the data
#==============================================================================

  save( oneRes, oneItRes, multiRes, file='../../data/res_flbeia.RData')
  
