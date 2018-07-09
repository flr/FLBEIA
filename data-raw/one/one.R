#-------------------------------------------------------------------------------  
# one dataset: 
# code to generate data in one.RData
# 
# Created: Agurtzane Urtizberea -  2016-02-09
# Changed: 2018-07-04 10:58:14 (ssanchez)
#------------------------------------------------------------------------------- 

# one.r - code to generate data in one.RData
# FLBEIA/data-row/one/one.r

# Copyright: AZTI, 2018
# Author: Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez (AZTI) (<flbeia@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


###############################################################################
# TITLE:        Test Case: 1 stock; 1 fleet; 1 season
#      
# NOTE #1:      Create FLBEIA input objects and run FLBEIA
###############################################################################

#-------------------------------------------------------------------------------
#   Section 1:       Load packages
#   Section 2:       Set working directory
#   Section 3:       Set input parameters related with time
#                         yrs <-  a vector with the next elements:
#                          first.yr      <- first year of historic data
#                          proj.yr       <- first year of projection
#                          last.yr       <- last year of projection
#   Section 4:       Set names, age, dimensions 
#                          fls           <- fleets' names                                    #vector
#                          stks          <- stocks' names                                    #vector
#                          fl1.mets      <- metiers' names in fleet 'fl1'                    #vector
#                          fl1.met1.stks  <- stocks' names in metier 'met1' and fleet 'fl1'    #vector
#                       all stocks:
#                           ni           <- number of iterations
#                           ns           <- number of seasons
#                       stock 'stk1':
#                           stk1.age.min  <- minimum age 
#                           stk1.age.max  <- maximum age 
#                           stk1.unit     <- number units 
#                       stock 'stk1'
#   Section 5:       biological data
#                       Historical data:  	
#                         stks.data <- a list with the name of stocks and 
#                                       in each stock all the data required:
#                           stk1_n.flq     <- Abundance at age     #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           stk1_m.flq     <- Natural mortality    #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           stk1_spwn.flq  <- Spawning             #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           stk1_fec.flq   <- Fecundity            #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           stk1_wt.flq    <- Weight               #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           stk1_range.min       <- minimum age
#                           stk1_range.max       <- maximum age
#                           stk1_range.plusgroup <- plusgroup
#                           stk1_range.minyear   <- minimum year
#                           stk1_range.minfbar   <- minimum age to take into account in the age 'f' calculation
#                           stk1_range.maxfbar   <- maximum age to take into account in the age 'f' calculation
#                       Projection parameters: in some variables, the projection is assumed the average of 
#                                               some historical years
#                           stk1_biol.proj.avg.yrs <- a vector with the years to calculate the average
#   Section 6:          Historical data per fleet:
#                         fleets.data <- a list with the name of the fleets and 
#                                       in each fleet all the data required:
#                         Fleet 'fl1':
#                           fl1.effort.flq       <- Effort                  #FLQuant[1,ny(hist),1,ns,1/ni]
#                           fl1.capacity.flq     <- Capacity                #FLQuant[1,ny(hist),1,ns,1/ni] (not required)
#                           fl1.fcost.flq        <- Fixed cost              #FLQuant[1,ny(hist),1,ns,1/ni] (not required)     
#                           fl1.crewshare.flq    <- Crewshare               #FLQuant[1,ny(hist),1,ns,1/ni] (not required)
#                         Data per fleet and metier:
#                          Fleet 'fl1' and metier 'met1'  
#                           fl1.met1.effshare.flq <- Effort share            #FLQuant[1,ny(hist),1,ns,1/ni] (not required)
#                         Data per fleet, metier and stock:
#                          Fleet 'fl1', metier 'met1' and stock 'stk1'  
#                           fl1.met1.stk1_landings.n.flq <- Landings at age                       #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           fl1.met1.stk1.discards.n     <- Discards at age                       #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                           fl1.met1.stk1_alpha.flq      <- alpha Cobb-Douglas parameter values   #FLQuant[na,ny,1/nu(stock),ns,1/ni]
#                           fl1.met1.stk1_beta.flq       <- beta Cobb-Douglas parameter values    #FLQuant[na,ny,1/nu(stock),ns,1/ni]
#                           fl1.met1.stk1_catch.q.flq    <- catch.q Cobb-Douglas parameter values #FLQuant[na,ny,1/nu(stock),ns,1/ni]
#                       Projection parameters: in some variables, the projection is assumed the average of some historical years
#                         Fleet 'fl1':
#                           fl1_proj.avg.yrs         <- a vector with the years to calculate the average,                           
#                                                          of the next variables; effort, capacity,crewshare
#                          Fleet 'fl1' and metier 'met1'  
#                           fl1.met1_proj.avg.yrs     <- a vector with the years to calculate the average,   
#                                                          of the next variables: vcost, effshare                       
#                          Fleet 'fl1', metier 'met1' and stock 'stk1'
#                           fl1.met1.stk1_proj.avg.yrs <- a vector with the years to calculate the average,
#                                                         of the next variables; price,landings.sel,discards.sel,
#                                                         landings.wt, discards.wt,alpha,beta,catch.q
#  Section 7:       SRs
#                   stks.data <- a list with the name of stocks and 
#                                       in each stock all the data required:
#                      stk1_sr.model        <- name of the SR model
#                      stk1_param.n         <- number of parameters in the SR model
#                      stk1_params.array    <- Parameters of the model      # array[param.n, ny, ns, ni]   
#                      stk1_params.name     <- Parameters' names            #vector
#                      stk1_rec.flq         <- Recruitment                  #FLQuant[1,ny(hist),1/nu(stock),ns,1/ni]
#                      stk1_ssb.flq         <- ssb                          #FLQuant[1,ny(hist),1/nu(stock),ns,1/ni]
#                      stk1_uncertainty.flq <- uncertainty of the model     #FLQuant[1,ny,1/nu(stock),ns,1/ni] (not required)
#                      stk1_proportion.flq  <- proportion of the model      #FLQuant[1,ny,1/nu(stock),ns,1/ni]
#                      stk1_prop.avg.yrs    <- a vector with the years to calculate the average of proportion
#                      stk1_timelag.matrix  <- timelag of year and season of spawning  #matrix(c(0,ni),2,ns, dimnames = list(c('year', 'season'),season = 1:ns))
#  Section 8:       advice:TAC/TAE/quota.share
#                     fleets object
#                     stks.data <- a list with the name of stocks and 
#                                       in each stock all the data required:
#                   
#                      stk1_advice.TAC.flq         <- TAC                          #FLQuant()        [1,ny,1,ns,1,ni]
#                      stk1_advice.quota.share.flq <- quota share between fleets   #FLQuant()        [1,ny,1,ns,1,ni]
#                      stk1_advice.avg.yrs         <- a vector with the years to calculate the average of (TAC/quota.share, if they are not defined)
#  Section 9:       indices
#  Section 10:       main.Ctrl
#                      main.ctrl           <- it's a list
#                      main.ctrl$sim.years <- a vector with the 'initial' and 'final' simulation year
#  Section 11:       biols.Ctrl
#                      growth.model     <- growth model 
#                      biols.ctrl       <- biols.ctrl object
#  Section 12:       fleets.Ctrl
#                      n.flts.stks      <- number of stks per fleet
#                      flts.stksnames   <- stocks name in each fleet
#                      effort.models    <- effort model
#                      effort.restr.fl1 <- stock restricting the effort in fl1 fleet 
#                      restriction.fl1  <- restricting the effort in fl1 fleet options: 'catch',...
#                      catch.models     <- catch model per fleet
#                      capital.models   <- capital model per fleet    
#                      flq              <- FLQuant with the dimenstions of stk1 stock. [1,ny,1,ns,1,ni] Required by fleets.ctrl #FLQuant()                           
#                      fleets.ctrl      <- fleets ctrl object  
#                      fleets.ctrl$fl1$stk1$discard.TAC.OS  <- logic value of discard.TAC.OS
#  Section 13:       advice.Ctrl
#                      HCR.models       <- HCR model per stock 
#                      ref.pts.stk1      <- reference points for 'stk1'
#                      advice.ctrl      <- advice ctrl object
#                      advice.ctrl[['stk1']][['sr']]            <- list() #defining the class of SR data for 'stk1' stock
#                      advice.ctrl[['stk1']][['sr']][['model']] <- character, model of SR data for 'stk1' stock
#                      advice.ctrl[['stk1']][['sr']][['years']] <- vector with y.rm value and num.yrs 
#  Section 14:       assess.Ctrl
#                      assess.models    <- assessment model per fleet
#                      assess.ctrl      <- assess.ctrl object
#                      assess.ctrl[['stk1']]$work_w_Iter     <- logical (TRUE if assessment model works with iters)
#                      assess.ctrl[['stk1']]$harvest.units   <- 'f' or 'hr'
#  Section 15:       obs.Ctrl
#                      stkObs.models    <- Observation model per stock
#                      flq.stk1          <- FLQuant with the dimenstions of stk1 stock. [1,ny,1,ns,1,ni] Required by obs.ctrl
#                      obs.ctrl         <- obs.ctrl object
#  Section 16:       BDs/covars/covars.ctrl/ NULL objects
#  Section 17:       FLBEIA input objects
#  Section 18:       Run FLBEIA
#  Section 19:       Results: Summary & Plot 
#  Section 20:       Change natural mortality rates in biols
#  Section 21:       Change stock-recruitment relationship
#  Section 22:       Change Fmsy in IcesHCR
#  Section 23:       Change fleets dynamic: Fixed Effort
#  Section 24:       Change HCR to annual TAC
#------------------------------------------------------------------------------#


#==============================================================================
# Section 1:            Install Packages
#==============================================================================

  library(FLCore) 
  library(FLAssess)
  library(FLash)
  library(FLFleet)
  library(FLXSA)
  library(FLBEIA) 


#==============================================================================
# Section 2:            Set working directory and Load the data
#==============================================================================

#wd <- ''
#setwd(wd)
  

#==============================================================================
# Section 3:            Set Simulation parameters related with time
#==============================================================================
  
  first.yr          <- 1990
  proj.yr           <- 2010 
  last.yr           <- 2025  
  yrs <- c(first.yr=first.yr,proj.yr=proj.yr,last.yr=last.yr)
  
#==============================================================================
#  Section 4:       Set names, age, dimensions 
#==============================================================================
  
  fls  <- c('fl1')

  stks <- c('stk1')

  fl1.mets      <- c('met1')
  fl1.met1.stks <- c('stk1')

  # all stocks the same
  
  ni <- 1
  it <- 1:ni
  ns <- 1
  
  # stock stk1
  stk1.age.min <- 1
  stk1.age.max <- 12
  stk1.unit    <- 1         
 
   
#==============================================================================
# Section 5:            biols
#==============================================================================
#  Data
#  stk1_n.flq, m, spwn, fec, wt
#==============================================================================

  # stock stk1
  stk1_n.flq     <- iter(as.FLQuant(read.csv(file = 'data/stk1_n.csv')),it)
  stk1_m.flq     <- iter(as.FLQuant(read.csv(file = 'data/stk1_m.csv')),it)
  stk1_spwn.flq  <- iter(as.FLQuant(read.csv(file = 'data/stk1_spwn.csv')),it)
  stk1_fec.flq   <- iter(as.FLQuant(read.csv(file = 'data/stk1_fec.csv')),it)
  stk1_wt.flq    <- iter(as.FLQuant(read.csv(file = 'data/stk1_wt.csv')),it)
  stk1_mat.flq   <- FLQuant(1,dimnames = list(age = stk1.age.min:stk1.age.max, year = first.yr:(proj.yr-1), unit = stk1.unit, season = 1:ns, iter = 1:ni))  

  stk1_range.min       <- 1
  stk1_range.max       <- 12
  stk1_range.plusgroup <- 12
  stk1_range.minyear   <- 1990
  stk1_range.minfbar   <- 4
  stk1_range.maxfbar   <- 9

  # Projection biols: weight,fecundity,mortality and spawning 
  
  #  we assume that the projection values of these variables are equal to 
  #  the average of the last 3 years of data.   
    
    stk1_biol.proj.avg.yrs <- c(2007:2009)
    
  #==============================================================================
  #              FLBEIA input object: biols
  #==============================================================================
  
    stks.data <- list(stk1=ls(pattern="^stk1")) 

    biols    <- create.biols.data(yrs,ns,ni,stks.data)
  
    # # plot biols
    # plotFLBiols(biols,pdfnm='s0')

    
#==============================================================================
# Section 6:            fleets
#==============================================================================
# Data per fleet
#    effort, crewshare, fcost, capacity
# Data per fleet and metier
#    effshare, vcost
# Data per fleet, metier and stock
#    landings.n, discards.n,landings.wt, discards.wt, landings, discards, landings.sel, discards.sel, price
#==============================================================================

  fl1_effort.flq        <- iter(as.FLQuant(read.csv(file = 'data/fl1_effort.csv')),it)
  fl1_capacity.flq      <- iter(as.FLQuant(read.csv(file = 'data/fl1_capacity.csv')),it)
  fl1_fcost.flq         <- FLQuant(417.05,dimnames = list(age = 'all', year = first.yr:(proj.yr-1), unit = stk1.unit, 
                                                           season = 1:ns, iter = 1:ni)) 
  fl1_crewshare.flq     <- FLQuant(0,dimnames = list(age = 'all', year = first.yr:(proj.yr-1), unit = stk1.unit, 
                                                          season = 1:ns, iter = 1:ni)) 

  fl1.met1_effshare.flq <- iter(as.FLQuant(read.csv(file = 'data/fl1.met1_effshare.csv')),it)
  fl1.met1_vcost.flq    <- FLQuant(0,dimnames = list(age = 'all', year = first.yr:(proj.yr-1), unit = stk1.unit, 
                                                                                               season = 1:ns, iter = 1:ni))  
  
  fl1.met1.stk1_landings.n.flq <- iter(as.FLQuant(read.csv(file = 'data/fl1.met1.stk1_landings.n.csv')),it)
  fl1.met1.stk1_discards.n.flq <- iter(as.FLQuant(read.csv(file = 'data/fl1.met1.stk1_discards.n.csv')),it)
  fl1.met1.stk1_price.flq <- FLQuant(0,dimnames = list(age = stk1.age.min:stk1.age.max, year = first.yr:(proj.yr-1), unit = stk1.unit, season = 1:ns, iter = 1:ni))  

 
  #         Projection
  #==============================================================================
  #         fleets: fl1
  #==============================================================================
  
	  fl1_proj.avg.yrs         <- c(2008:2009)
  	fl1.met1_proj.avg.yrs      <- c(2008:2009)   
  	fl1.met1.stk1_proj.avg.yrs <- c(2006:2008)   

  #==============================================================================
  #              FLBEIA input object: fleets
  #==============================================================================
  
  	fls.data <- list(fl1=ls(pattern="^fl1")) 
  	
  	fleets   <- create.fleets.data(yrs,ns,ni,fls.data,stks.data)
  	
  	fleets[["fl1"]]@metiers[["met1"]]@catches[["stk1"]]@price[] <- rnorm(length(stk1.age.min:stk1.age.max)*length(first.yr:last.yr)*ns*ni, 2500,250)
  	fleets[["fl1"]]@crewshare[] <- rnorm(length(first.yr:last.yr)*ns*ni, .25,.03)
  	fleets[["fl1"]]@metiers[["met1"]]@vcost[] <-rnorm(length(first.yr:last.yr)*ns*ni, 1000,100) 
  	
  	# # plot fleets
  	# plotFLFleets(fleets,pdfnm='s0')
  	

#==============================================================================
#  Section 7:       SRs
#==============================================================================

  stk1_sr.model        <- 'bevholtAR1'
  stk1_params.n        <- 3
  stk1_params.array    <- xtabs2(data~param+year+season+iter, 
                                data=read.csv(file = 'data/stk1_params.csv'),
                                exclude=NULL,na.action=na.pass)[,,,it,drop=F]         
  stk1_params.name     <- c('a','b','c') 
  stk1_rec.flq         <- iter(as.FLQuant(read.csv(file = 'data/stk1_rec.csv')),it)
  stk1_ssb.flq         <- iter(as.FLQuant(read.csv(file = 'data/stk1_ssb.csv')),it)
  stk1_uncertainty.flq <- iter(as.FLQuant(read.csv(file = 'data/stk1_uncertainty.csv')),it)
  stk1_proportion.flq  <- iter(as.FLQuant(read.csv(file = 'data/stk1_proportion.csv')),it)
  stk1_prop.avg.yrs    <- ac(2006:2008)
  stk1_timelag.matrix  <- matrix(c(1,1),nrow=2,ncol=1, dimnames = list(c('year', 'season'),'all'))

  #==============================================================================
  #              FLBEIA input object: SRs
  #==============================================================================
  
    stks.data <- list(stk1=ls(pattern="^stk1")) 
    
    SRs      <- create.SRs.data(yrs,ns,ni,stks.data)
  
  
#==============================================================================
#  Section 8:       advice:TAC/TAE/quota.share
#==============================================================================

  stk1_advice.TAC.flq         <- iter(as.FLQuant(read.csv(file = 'data/stk1_advice.tac.csv')),it)
  stk1_advice.quota.share.flq <- iter(as.FLQuant(read.csv(file = 'data/stk1_advice.quota.share.csv')),it)
  stk1_advice.avg.yrs         <- c(2006:2008)

  #==============================================================================
  #              FLBEIA input object: advice
  #==============================================================================
  
    stks.data <- list(stk1=ls(pattern="^stk1")) 
  
    advice   <- create.advice.data(yrs,ns,ni,stks.data,fleets)
    
  
#==============================================================================
#  Section 9:       indices
#==============================================================================

  indices <- NULL
  
  #.................based on biomass.......................
  
  flq   <- quantSums(biols[[1]]@n*biols[[1]]@wt)
  unc   <- id <- q <- flq
  unc[] <- rlnorm( prod(dim(flq)), 0, 0.3)
  q[]   <- rep( runif(dim(flq)[1], 1e-05/2, 1e-05*5), dim(flq)[2])
  id    <- flq*unc*q
  
  stk1_indices <- "idBio"
  stk1_idBio_index.flq  <- id
  stk1_idBio_index.q.flq <- q
  stk1_idBio_index.var.flq <- unc
  stk1_idBio_range.startf <- 0.12
  stk1_idBio_range.endf <- 1-0.12
  stk1_idBio_range.min <- 0
  stk1_idBio_range.max <- 0
  #  YFT_cpue_range.plusgroup <- 0
  stk1_idBio_range.minyear <- first.yr
  stk1_idBio_range.maxyear <- proj.yr-1
  
  stk1_idBio_type <- "FLIndexBiomass"
  stks.data <- list(stk1=ls(pattern="^stk1")) 
  indices.bio <- create.indices.data(yrs,ns,ni,stks.data)
  
  #.................based on age.......................
  library(gtools)
  flq   <- biols[["stk1"]]@n
  unc   <- id <- q <- flq
  unc[] <- rlnorm(prod(dim(flq)), 0,0.3)
  q[]   <- rep(runif(dim(flq)[1], 1e-05/2, 1e-05*5), dim(flq)[2])
  id    <- biols[["stk1"]]@n*unc*q

  stk1_indices <- "idAge"
  stk1_idAge_index.flq  <- id
  stk1_idAge_index.q.flq <- q
  stk1_idAge_index.var.flq <- unc
  stk1_idAge_range.startf <- 0.12
  stk1_idAge_range.endf <- 1-0.12
  stk1_idAge_range.min <- stk1.age.min
  stk1_idAge_range.max <- stk1.age.max
  #  YFT_cpue_range.plusgroup <- 0
  stk1_idAge_range.minyear <- first.yr
  stk1_idAge_range.maxyear <- proj.yr-1
  
  stk1_idAge_type <- "FLIndex"
  stks.data <- list(stk1=ls(pattern="^stk1")) 
  indices.age <- create.indices.data(yrs,ns,ni,stks.data)
  

#==============================================================================
#  Section 10:       main.ctrl
#==============================================================================

  main.ctrl           <- list()
  main.ctrl$sim.years <- c(initial = proj.yr, final = last.yr)

  
#==============================================================================
#  Section 11:       biol.ctrl
#==============================================================================

  growth.model     <- c('ASPG')
  biols.ctrl       <- create.biols.ctrl (stksnames=stks,growth.model=growth.model)

  
#==============================================================================
#  Section 12:       fleets.ctrl
#==============================================================================

  n.fls.stks       <- 1
  fls.stksnames    <- 'stk1'
  effort.models    <- 'SMFB'
  effort.restr.fl1 <- 'stk1'
  restriction.fl1  <- 'catch'
  catch.models     <- 'CobbDouglasAge'
  capital.models   <- 'fixedCapital'       
  flq.stk1<- FLQuant(dimnames = list(age = 'all', year = first.yr:last.yr, unit = stk1.unit, 
                                     season = 1:ns, iter = 1:ni)) 
  fleets.ctrl      <- create.fleets.ctrl(fls=fls,n.fls.stks=n.fls.stks,fls.stksnames=fls.stksnames,
                                         effort.models= effort.models, catch.models=catch.models,
                                         capital.models=capital.models, flq=flq.stk1,
                                         effort.restr.fl1 = effort.restr.fl1, restriction.fl1 = restriction.fl1)

  fleets.ctrl$fl1$stk1$discard.TAC.OS  <- FALSE
  fleets.ctrl$fl1$restriction <- "landings"


#==============================================================================
#  Section 13:       advice.ctrl
#==============================================================================

  HCR.models   <- c('IcesHCR') 
  ref.pts.stk1 <- matrix(rep(c(548.6296271, 768.0814779, 0.1057783),3), 3, ni, dimnames = list(c('Blim', 'Btrigger','Fmsy'), 1:ni))
  advice.ctrl  <- create.advice.ctrl(stksnames = stks, HCR.models =  HCR.models, 
                                     ref.pts.stk1 = ref.pts.stk1,first.yr=first.yr,last.yr=last.yr)
    
  advice.ctrl[['stk1']][['sr']]            <- list()
  advice.ctrl[['stk1']][['sr']][['model']] <- 'geomean'
  advice.ctrl[['stk1']][['sr']][['years']] <- c(y.rm = 2, num.years = 10)
  advice.ctrl$stk1$AdvCatch <- rep(TRUE,length(first.yr:last.yr))   #TRUE advice in catches, FALSE advice in landings
  names(advice.ctrl$stk1$AdvCatch) <- as.character((first.yr:last.yr))
  

#==============================================================================
#  Section 14:       assess.ctrl
#==============================================================================
  
  assess.models    <- 'NoAssessment'
  assess.ctrl <- create.assess.ctrl(stksnames = stks, assess.models = assess.models)
  assess.ctrl[['stk1']]$work_w_Iter   <- TRUE
  
  #==============================================================================
  #source("sca.wrapper.R")
  
  assess.ctrl.age <- assess.ctrl
  assess.ctrl.age$stk1$assess.model <- "sca2flbeia"
  assess.ctrl.age[["stk1"]]$harvest.units <- "f"
  assess.ctrl.age[["stk1"]]$control$test <- TRUE
  
  #==============================================================================
  
  assess.ctrl.bio <- assess.ctrl.age
  assess.ctrl.bio[["stk1"]]$assess.model <- "spict2flbeia"
  

#==============================================================================
#  Section 15:       obs.ctrl
#==============================================================================
  
  stkObs.models <- "perfectObs"
  flq.stk1 <- FLQuant(dimnames = list(age = 'all', year = first.yr:last.yr, unit = stk1.unit, 
                                      season = 1:ns, iter = 1:ni)) 
  
  obs.ctrl <- create.obs.ctrl(stksnames = stks,  stkObs.models = stkObs.models,flq.stk1 = flq.stk1)
  
  #............................
  
  #==========================OBSERVATION: AGE2BIODAT====================================================
  
  stkObs.models <- 'age2bioDat'
  flq.stk1 <- FLQuant(dimnames = list(age = 'all', year = first.yr:last.yr, unit = stk1.unit, 
                                      season = ns, iter = 1:ni)) 
  
  obs.ctrl.bio <- create.obs.ctrl(stksnames = stks,  stkObs.models = stkObs.models, flq.stk1=flq.stk1)
  
  obs.ctrl.bio[['stk1']][['indObs']] <- vector('list', 1)
  names(obs.ctrl.bio[['stk1']][['indObs']]) <- c("idBio")
  
  obs.ctrl.bio[['stk1']][['indObs']][['idBio']] <- list()
  obs.ctrl.bio[['stk1']][['indObs']][['idBio']][['indObs.model']]  <- 'bioInd'
  obs.ctrl.bio$stk1$stkObs$TAC.ovrsht <- array(1,dim=c(length(stks),length(first.yr:last.yr)),dimnames=list(c(stks),ac(first.yr:last.yr)))
  obs.ctrl.bio$stk1$stkObs$land.bio.error <- array(1,dim=c(length(stks),length(first.yr:last.yr)),dimnames=list(c(stks),ac(first.yr:last.yr)))
  obs.ctrl.bio$stk1$stkObs$disc.bio.error <- array(1,dim=c(length(stks),length(first.yr:last.yr)),dimnames=list(c(stks),ac(first.yr:last.yr)))

  
  #=========================OBSERVATION: AGE2AGEDAT================================================
  
    flq <- biols[["stk1"]]@n
    na  <- stk1.age.max
    ny  <- length(first.yr:last.yr)
    
    ages.error <- array(0,dim = c(na, na, ny,ni))
    for(a in 1:12){
      for(i in 1:ni){
        for(y in 1:ny){
          if(a == 1)  ages.error[1,,y,i]  <- rdirichlet(1, c(0.85,0.1, 0.05, rep(0,9)))
          if(a == 2)  ages.error[2,,y,i]  <- rdirichlet(1, c(0.1,0.75,0.1, 0.05, rep(0,8)))
          if(a == 3)  ages.error[3,,y,i]  <- rdirichlet(1, c(0.05,0.1,0.7,0.1, 0.05, rep(0,7)))
          if(a == 4)  ages.error[4,,y,i]  <- rdirichlet(1, c(rep(0,1), 0.05,0.1,0.7,0.1,0.05, rep(0,6)))
          if(a == 5)  ages.error[5,,y,i]  <- rdirichlet(1, c(rep(0,2), 0.05,0.1,0.7,0.1,0.05, rep(0,5)))
          if(a == 6)  ages.error[6,,y,i]  <- rdirichlet(1, c(rep(0,3), 0.05,0.1,0.7,0.1,0.05, rep(0,4)))
          if(a == 7)  ages.error[7,,y,i]  <- rdirichlet(1, c(rep(0,4), 0.05,0.1,0.7,0.1,0.05, rep(0,3)))
          if(a == 8)  ages.error[8,,y,i]  <- rdirichlet(1, c(rep(0,5), 0.05,0.1,0.7,0.1,0.05, rep(0,2)))
          if(a == 9)  ages.error[9,,y,i]  <- rdirichlet(1, c(rep(0,6), 0.05,0.1,0.7,0.1,0.05, rep(0,1)))
          if(a == 10) ages.error[10,,y,i] <- rdirichlet(1, c(rep(0,7), 0.05,0.1,0.7,0.1,0.05))
          if(a == 11) ages.error[11,,y,i] <- rdirichlet(1, c(rep(0,8), 0.05,0.1,0.75,0.1))
          if(a == 12) ages.error[12,,y,i] <- rdirichlet(1, c(rep(0,9), 0.05,0.1,0.85))
        }                                                                                      
      }} 
    
    nmort.error      <- fec.error       <- land.wgt.error  <- stk.nage.error <- stk.wgt.error <-
      disc.wgt.error <- land.nage.error <- disc.nage.error <- flq
     
    TAC.ovrsht <- flq[1,]
    dimnames(TAC.ovrsht)[[1]] <- 'all' 
  
  #=========================OBSERVATION: AGE2AGEPOP================================================
  
    # age2agePop  oneObsCStk
  
    stkObs.models <- 'age2ageDat'
    flq.stk1 <- FLQuant(dimnames = list(age = 'all', year = first.yr:last.yr, unit = stk1.unit, 
                                        season = ns, iter = 1:ni)) 
    
    obs.ctrl.age <- create.obs.ctrl(stksnames = stks, stkObs.models = stkObs.models, flq.stk1=flq.stk1)
    obs.ctrl.age[['stk1']][['indObs']] <- vector('list', 1)
    names(obs.ctrl.age[['stk1']][['indObs']]) <- c("idAge")
    
    obs.ctrl.age[['stk1']][['indObs']][['idAge']] <- list()
    obs.ctrl.age[['stk1']][['indObs']][['idAge']][['indObs.model']]  <- 'ageInd'
    obs.ctrl.age$stk1$stkObs$TAC.ovrsht <- TAC.ovrsht
    # obs.ctrl.age$stk1$stkObs$TAC.ovrsht <- array(1,dim=c(length(stks),length(first.yr:last.yr)),dimnames=list(c(stks),ac(first.yr:last.yr)))
    # obs.ctrl.age$stk1$stkObs$land.bio.error <- array(1,dim=c(length(stks),length(first.yr:last.yr)),dimnames=list(c(stks),ac(first.yr:last.yr)))
    # obs.ctrl.age$stk1$stkObs$disc.bio.error <- array(1,dim=c(length(stks),length(first.yr:last.yr)),dimnames=list(c(stks),ac(first.yr:last.yr)))
    obs.ctrl.age[['stk1']][['stkObs']][['ages.error']] <- ages.error
    
    slts <- c('nmort.error', 'fec.error', 'land.wgt.error', 'stk.nage.error', 'stk.wgt.error',
              'disc.wgt.error', 'land.nage.error', 'disc.nage.error')
    for(sl in slts){
      obs.ctrl.age[['stk1']][['stkObs']][[sl]]   <- get(sl)
      obs.ctrl.age[['stk1']][['stkObs']][[sl]][] <- rnorm(prod(dim(flq)),1,.1)
    }
  
  
#==============================================================================
#  Section 17:       covars
#==============================================================================

  cv_mean_value <- c( FuelCost = 46, CapitalCost = 4519.06, Salaries = 0, InvestShare = 0.2, NumbVessels = 228.33, 
                      MaxDays = 228, w1 = 0.03, w2 = 0.03, EmploymentPerVessel = 2)

  covars <- vector("list",9)
  names(covars) <- c("FuelCost","CapitalCost","Salaries", "InvestShare","NumbVessels","MaxDays",
                     "w1","w2","EmploymentPerVessel")
  
  flq <- FLQuant( rnorm(length(fls)*length(first.yr:last.yr)*ns*ni, 1000,100), 
                  dimnames = list(fleets = fls, year = first.yr:last.yr, unit = stk1.unit, season = 1:ns, iter = 1:ni)) 
  
  for(cv in names(covars)){ 
    covars[[cv]] <- flq
    dimnames(covars[[cv]])[[1]] <- c('fl1')
    # The mean value is cv_mean_value[cv] and the CV is equal to 10%.
    covars[[cv]][] <- rnorm(prod(dim(flq)), cv_mean_value[cv], cv_mean_value[cv]*0.1)
  }
  
  
#==============================================================================
#  Section 18:       covars.ctrl
#==============================================================================
  
  covars.ctrl <- vector("list",9)
  names(covars.ctrl) <- c("FuelCost","CapitalCost","Salaries", "InvestShare","NumbVessels","MaxDays", 
                          "w1","w2","EmploymentPerVessel")
  
  for(cv in names(covars)){ 
    covars.ctrl[[cv]] <- list()
    covars.ctrl[[cv]][['process.model']]  <- 'fixedCovar'
  }
  
  
#==============================================================================
#  Section 16:       BDs
#==============================================================================
  
  BDs <- NULL
  
  
#==============================================================================
#  Section 17:       Rename
#==============================================================================
  
  oneAdv        <- advice
  oneAdvC       <- advice.ctrl
  oneAssC       <- assess.ctrl
  oneBio        <- biols
  oneBioC       <- biols.ctrl
  oneCv         <- covars
  oneCvC        <- covars.ctrl
  oneFl         <- fleets
  oneFlC        <- fleets.ctrl
  oneMainC      <- main.ctrl
  oneObsC       <- obs.ctrl
  oneSR         <- SRs
  oneIndBio     <- indices.bio
  oneIndAge     <- indices.age
  oneObsC       <- obs.ctrl
  oneObsCIndBio <- obs.ctrl.bio
  oneObsCIndAge <- obs.ctrl.age


#==============================================================================
#  Section 17:       FLBEIA input objects
#==============================================================================
  
  save( oneAdv, oneAdvC, oneAssC, oneBio, oneBioC, oneCv, oneCvC, 
        oneFl, oneFlC, oneIndAge, oneIndBio, 
        oneMainC, oneObsC, oneObsCIndAge, oneObsCIndBio, oneSR, 
        file ="../../data/one.RData")
  
  
# #==============================================================================
# #  Section 18:       Run FLBEIA
# #==============================================================================
# 
#   s0 <- FLBEIA(biols = oneBio,    # FLBiols object with one FLBiol element for stk1.
#                SRs = oneSR,       # A list with one FLSRSim object for stk1.
#                BDs = NULL,        # No Biomass Dynamic populations in this case.
#                fleets = oneFl,    # FLFleets object with on fleet.
#                covars = oneCv,    # covars not used
#                indices = NULL,    # indices not used 
#                advice = oneAdv,   # A list with two elements 'TAC' and 'quota.share'
#                main.ctrl = oneMainC,   # A list with one element to define the start and end of the simulation.
#                biols.ctrl = oneBioC,   # A list with one element to select the model to simulate the stock dynamics.
#                fleets.ctrl = oneFlC,   # A list with several elements to select fleet dynamic models and store additional parameters.
#                covars.ctrl = oneCvC,   # covars control not used 
#                obs.ctrl = oneObsC,     # A list with one element to define how the stock observed ("PerfectObs").
#                assess.ctrl = oneAssC,  # A list with one element to define how the stock assessment model used ("NoAssessment").
#                advice.ctrl = oneAdvC)  # A list with one element to define how the TAC advice is obtained ("IcesHCR").
# 
#   oneRes <- s0
#   save(oneRes,file='../data/oneRes.RData')
#   
#   
# #==============================================================================
# #  Section 19:       Results: Summary & Plot 
# #==============================================================================
# 
# # Names of the object returned by FLBEIA
# names(s0)
# 
# # The default plot for FLBiol defined in FLCore
# plot(s0$biols[[1]])
# 
# # Create summary data frames (biological, economic, and catch)
# 
# proj.yr     <- 2013 
# s0_sum      <- bioSum(s0)
# s0_flt      <- fltSum(s0)
# s0_fltStk   <- fltStkSum(s0)
# 
# 
# # Create several plots and save them in the working directory using 'pdf' format and 
# # 's0' suffix in the name.
# 
# plotFLBiols(s0$biols, pdfnm='oneRes')
# plotFLFleets(s0$fleets,pdfnm='oneRes')
# plotfltStkSum(s0, pdfnm='oneRes') 
# plotEco(s0$fleets,pdfnm='oneRes')

