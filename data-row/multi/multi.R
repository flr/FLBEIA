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
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Test Case: 2stock and 2 fleets
# NOTE #1:     
# NOTE #2:      Create FLBEIA input objects and run FLBEIA
###############################################################################

#-------------------------------------------------------------------------------
#   Section 1:       Charge libraries
#   Section 2:       Set directory
#   Section 3:       Load the data
#   Section 4:       Set Simulation-time parameters
#                         yrs <-  a vector with the next elements:
#                          first.yr      <- first year of historic data
#                          proj.yr       <- first year of projection
#                          last.yr       <- last year of projection
#   Section 5:       Set names, age, dimensions 
#                          fls           <- fleets' names                                    #vector
#                          stks          <- stocks' names                                    #vector
#                          fl1.mets      <- metiers' names in fleet 'fl1'                    #vector
#                          fl1.met1.stks <- stocks' names in metier 'met1' and fleet 'fl1'   #vector
#                       all stocks:
#                           ni           <- number of iterations
#                           ns           <- number of seasons
#                       stock 'stk1':
#                           stk1.age.min  <- minimum age 
#                           stk1.age.max  <- maximum age
#                           stk1.unit     <- number units
#                       stock 'stk2':
#   Section 6:       biological data
#                       Historical data:
#                         Stock 'stk1':
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
#                         stock 'stk2': 
#                       Projection parameters: in some variables, the projection is assumed the average of some historical years
#                         stock 'stk1': 
#                           stk1_biol.proj.avg.yrs <- a vector with the years to calculate the average
#                         stock 'stk2': 
#   Section 7:          Historical data per fleet:
#                         fleets.data <- a list with the name of the fleets and 
#                                       in each fleet all the data required:
#                         Fleet 'fl1' 
#                           fl1.effort.flq        <- Effort          #FLQuant[1,ny(hist),1,ns,1/ni]
#                           fl1.capacity.flq      <- Capacity        #FLQuant[1,ny(hist),1,ns,1/ni](not required)
#                         Historical data per fleet and metier:
#                           Fleet 'fl1' and metier 'met1' 
#                             fl1.met1.effshare.flq <-Effort share                  #FLQuant[1,ny(hist),1,ns,1/ni]
#                             fl1.met1.vcost.flq    <- Variable cost                #FLQuant[1,ny(hist),1,ns,1/ni](not required)
#                           Historical data per fleet,metier and stock:
#                             Fleet 'fl1', metier 'met1' and stock 'stk1':
#                               fl1.met1.stk1_landings.n.flq    <- Landings at age of the stock      #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                               fl1.met1.stk1_discards.n.flq    <- Discards at age of the stock      #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni]
#                               fl1.met1.stk1_price.flq         <- Price of the stock                #FLQuant[na,ny(hist),1/nu(stock),ns,1/ni] (not required)
#                               fl1.met1.stk1_alpha.flq     <- alpha Cobb-Douglas parameter values   #FLQuant[na,ny,1/nu(stock),ns,1/ni] (not required)
#                               fl1.met1.stk1_beta.flq      <- beta Cobb-Douglas parameter values    #FLQuant[na,ny,1/nu(stock),ns,1/ni] (not required)
#                               fl1.met1.stk1_catch.q.flq   <- catch.q Cobb-Douglas parameter values #FLQuant[na,ny,1/nu(stock),ns,1/ni] (not required)
#                       Projection parameters: in some variables, the projection is assumed the average of some historical years
#                         fleet 'fl1'
#                           fl1_proj.avg.yrs           <-  a vector with the years to calculate the average,            
#                                                          of the next variables; effort, capacity,crewshare      
#                           fleet 'fl1' and metier 'met1':
#                           fl1.met1_proj.avg.yrs      <- a vector with the years to calculate the average,   
#                                                          of the next variables: vcost, effshare                       
#                           fl1.met1.stk1_proj.avg.yrs <- a vector with the years to calculate the average,
#                                                         of the next variables; price,landings.sel,discards.sel,
#                                                         landings.wt, discards.wt,alpha,beta,catch.q
#                         fleet 'fl2'
#  Section 8:       SRs
#                   stks.data <- a list with the name of stocks and 
#                                       in each stock all the data required:
#                      stk1_sr.model        <- name of the SR model
#                      stk1_param.n         <- number of parameters in the SR model
#                      stk1_params.array    <- Parameters of the model      # array[param.n, ny, ns, ni]   
#                      stk1_params.name     <- Parameters' names
#                      stk1_rec.flq         <- Recruitment                  #FLQuant[1,ny(hist),1/nu(stock),ns,1/ni]
#                      stk1_ssb.flq         <- ssb                          #FLQuant[1,ny(hist),1/nu(stock),ns,1/ni]
#                      stk1_uncertainty.flq <- uncertainty of the model     #FLQuant[na,ny,1/nu(stock),ns,1/ni] (not required)
#                      stk1_proportion.flq  <- proportion of the model      #FLQuant[na,ny,1/nu(stock),ns,1/ni]
#                      stk1_prop.avg.yrs    <- a vector with the years to calculate the average of proportion
#                      stk1_timelag.matrix  <- timelag of year and season of spawning  #matrix(c(0,1),2,ns, dimnames = list(c('year', 'season'),season = 1:ns)) 
#  Section 9:       BDs
#                      stk2_bd.model        <- Biomass dynamic model
#                      stk2_param.n         <- Number of parameters in the model
#                      stk2_params.array    <- Parameters of the model          # array[param.n, ny, ns, ni]   
#                      stk2_params.name     <- Parameters' names
#                      stk2_biomass.flq     <- Biomass                          #FLQuant[1,ny(hist),1/nu(stock),ns,1/ni]
#                      stk2_catch.flq       <- Catch                            #FLQuant[1,ny(hist),1/nu(stock),ns,1/ni]
#                      stk2_uncertainty.flq <- Uncertainty of the model         #FLQuant[1,ny,1/nu(stock),ns,1/ni]   (not required)
#  Section 10:       advice:TAC/TAE/quota.share   
#                     fleets object
#                     stks.data <- a list with the name of stocks and 
#                                       in each stock all the data required:
#                      stk1_advice.TAC.flq  <- TAC                #FLQuant[1,ny,1,1,1/ni]
#                      Projection:
#                      stk1_advice.avg.yrs  <-a vector with the years to calculate the average of (TAC/quota.share, if they are not defined)
#  Section 11:       main.ctrl
#                      main.ctrl           <- it's a list
#                      main.ctrl$sim.years <- a vector with the 'initial' and 'final' simulation year
#  Section 12:       biols.ctrl
#                      growth.model     <- growth model 
#                      biols.ctrl       <- biols.ctrl object
#  Section 13:       fleets.ctrl
#                      n.flts.stks      <- number of stks per fleet
#                      flts.stksnames   <- stocks' names in each fleet
#                      effort.models    <- effort model per fleet
#                      effort.restr.fl2 <- stock restricting the effort in fl2 fleet 
#                      restriction.fl2  <- restricting the effort in fl2 fleet options: 'catch',...
#                      catch.models     <- catch  model per fleet
#                      capital.models   <- capital model per fleet
#                      price.models     <- elastic price  
#                      flq     <- FLQuant with the dimenstions of the stock. Required by fleets.ctrl                 
#                      fleets.ctrl      <- fleets ctrl object  
#  Section 14:       advice.ctrl
#                      HCR.models       <- HCR model per stock 
#                      Ftarget.stk2     <- 'stk2' stock's f target
#                      ref.pts.stk1     <- 'stk1' stock's reference points
#                      advice.ctrl      <- advice ctrl object
#  Section 15:       assess.ctrl
#                      assess.models    <- assessment model per fleet
#                      assess.ctrl      <- assess.ctrl object
#  Section 16:       obs.ctrl
#                      stkObs.models    <- Observation model per stock
#                      obs.ctrl         <- obs.ctrl object
#  Section 17:       covars.ctrl
#                      covars.ctrl      <- NULL
#  Section 18:       FLBEIA input objects
#  Section 19:       Run FLBEIA
#  Section 20:       Plot
#------------------------------------------------------------------------------#

rm(list=ls())
# library(devtools)
# install.packages(c("plyr", "ggplot2", "nloptr", "mvtnorm", "triangle"))
# install.packages("flr/FLCore")
# install.packages("flr/FLFleet")
# install.packages(c("FLAssess","FLash","FLXSA"), repos="http://flr-project.org/R")
# install_github("flr/FLBEIA")


#==============================================================================
# Section 1:            Charge libraries
#==============================================================================

  library(FLCore) 
  library(FLAssess)
  library(FLash)
  library(FLFleet)
  library(FLXSA)
  library(FLBEIA) 

#==============================================================================
# Section 2:            Set directory and Load the data
#==============================================================================

#wd <- ''
#setwd(wd)


#==============================================================================
# Section 3:            Load the data
#==============================================================================

  #load('2stock2fleets4seasons/data/input_data.RData')


#==============================================================================
# Section 4:            Set Simulation parameters related with time
#==============================================================================
  
  first.yr <- 1990
  proj.yr  <- 2009 
  last.yr  <- 2025  
  yrs <- c(first.yr=first.yr,proj.yr=proj.yr,last.yr=last.yr)
  
  
#==============================================================================
#  Section 5:       Set names, age, dimensions 
#==============================================================================
  
  fls   <- c('fl1','fl2')

  stks <- c('stk1','stk2')

  fl1.mets      <- c('met1','met2')
  fl1.met1.stks <- c('stk1','stk2')
  fl1.met2.stks <- c('stk1','stk2')

  fl2.mets      <- c('met1','met2')
  fl2.met1.stks <- c('stk1','stk2')
  fl2.met2.stks <- c('stk1','stk2')

  # all stocks the same
  ni <- 1
  ns <- 4
  
  # stock stk1
  stk1.age.min <- 0
  stk1.age.max <- 15
  stk1.unit    <- 4         
  
  # stock stk2
  stk2.age.min <- NA
  stk2.age.max <- NA
  stk2.unit    <- 1         

  
#==============================================================================
# Section 6:            biols
#==============================================================================
#  Historical data
#  stk1_n.flq, m, spwn, fec, wt
#==============================================================================

  #stock stk1
  stk1_n.flq     <- as.FLQuant(read.csv(file = 'data/stk1_n.csv'))
  stk1_m.flq     <- as.FLQuant(read.csv(file = 'data/stk1_m.csv'))
  stk1_spwn.flq  <- as.FLQuant(read.csv(file = 'data/stk1_spwn.csv'))
  stk1_mat.flq   <- as.FLQuant(read.csv(file = 'data/stk1_mat.csv'))
  stk1_fec.flq   <- stk1_mat.flq
  stk1_fec.flq[] <- 1
  stk1_wt.flq    <- as.FLQuant(read.csv(file = 'data/stk1_wt.csv'))
  
  stk1_range.min       <- 0
  stk1_range.max       <- 15
  stk1_range.plusgroup <- 15
  stk1_range.minyear   <- 1990
  stk1_range.minfbar   <- 1
  stk1_range.maxfbar   <- 5

  # stock stk2
  stk2_n.flq     <- as.FLQuant(read.csv(file = 'data/stk2_n.csv'))
  stk2_m.flq     <- as.FLQuant(read.csv(file = 'data/stk2_m.csv'))
  stk2_spwn.flq  <- as.FLQuant(read.csv(file = 'data/stk2_spwn.csv'))
  stk2_mat.flq   <- as.FLQuant(read.csv(file = 'data/stk2_mat.csv'))
  stk2_fec.flq   <- stk2_mat.flq
  stk2_fec.flq[] <- 1
  stk2_wt.flq    <- as.FLQuant(read.csv(file = 'data/stk2_wt.csv'))

  stk2_range.min       <- NA
  stk2_range.max       <- NA
  stk2_range.plusgroup <- NA
  stk2_range.minyear   <- 1990
  stk2_range.minfbar   <- 1
  stk2_range.maxfbar   <- 1

  #               Projection: 
  #  we assume that the projection values of some variables are equal to 
  #  the average of some historical years:weight,fecundity,mortality and spawning    
  
  stk1_biol.proj.avg.yrs <- c(2006:2008)
  stk2_biol.proj.avg.yrs <- c(2006:2008)

    
  #==============================================================================
  #              FLBEIA input object: biols
  #==============================================================================
  
    stks.data <- list(stk1=ls(pattern="^stk1"),stk2=ls(pattern="^stk2")) 
  
    biols <- create.biols.data(yrs,ns,ni,stks.data)
  
    # # plot biols
    # plotFLBiols(biols,pdfnm='s0')


#==============================================================================
# Section 7:            fleets
#==============================================================================
# Data per fleet
#    effort, crewshare, fcost, capacity
# Data per fleet and metier
#    effshare, vcost
# Data per fleet, metier and stock
#    landings.n, discards.n,landings.wt, discards.wt, landings, discards, landings.sel, discards.sel, price
#==============================================================================
        
  fl1_effort.flq        <- as.FLQuant(read.csv(file = 'data/fl1_effort.csv'))
  fl1_capacity.flq      <- as.FLQuant(read.csv(file = 'data/fl1_capacity.csv'))
  fl1.met1_effshare.flq <- as.FLQuant(read.csv(file = 'data/fl1.met1_effshare.csv'))
  fl1.met2_effshare.flq <- as.FLQuant(read.csv(file = 'data/fl1.met2_effshare.csv'))
  
  fl1.met1.stk1_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl1.met1.stk1_landings.n.csv'))
  fl1.met1.stk1_discards.n.flq <- as.FLQuant(read.csv(file = 'data/fl1.met1.stk1_discards.n.csv'))
  fl1.met1.stk2_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl1.met1.stk2_landings.n.csv'))
  fl1.met1.stk2_discards.n.flq <- as.FLQuant(read.csv(file = 'data/fl1.met1.stk2_discards.n.csv'))
  
  fl1.met2.stk1_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl1.met2.stk1_landings.n.csv'))
  fl1.met2.stk2_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl1.met2.stk2_landings.n.csv'))
  
  fl2_effort.flq   <- as.FLQuant(read.csv(file = 'data/fl2_effort.csv'))
  fl2_capacity.flq <- as.FLQuant(read.csv(file = 'data/fl2_capacity.csv'))
  fl2_fcost.flq    <- as.FLQuant(read.csv(file = 'data/fl2_fcost.csv'))
  
  fl2.met1_effshare.flq <- as.FLQuant(read.csv(file = 'data/fl2.met1_effshare.csv'))
  fl2.met2_effshare.flq <- as.FLQuant(read.csv(file = 'data/fl2.met2_effshare.csv')) 
  fl2.met1_vcost.flq    <- as.FLQuant(read.csv(file = 'data/fl2.met1_vcost.csv'))
  fl2.met2_vcost.flq    <- as.FLQuant(read.csv(file = 'data/fl2.met2_vcost.csv'))
  
  fl2.met1.stk1_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met1.stk1_landings.n.csv'))
  fl2.met2.stk1_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met2.stk1_landings.n.csv'))
  fl2.met1.stk2_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met1.stk2_landings.n.csv'))
  fl2.met2.stk2_landings.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met2.stk2_landings.n.csv'))
  fl2.met1.stk1_discards.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met1.stk1_discards.n.csv'))
  fl2.met2.stk1_discards.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met2.stk1_discards.n.csv'))
  fl2.met1.stk2_discards.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met1.stk2_discards.n.csv'))
  fl2.met2.stk2_discards.n.flq <- as.FLQuant(read.csv(file = 'data/fl2.met2.stk2_discards.n.csv'))
  
  # Cobb Douglas parameters estimates by default in the function create.fleets.data
  # when the stock is age structured
  #       alpha=1, beta=1,q=Catches/N*E
  #
  # Cobb Douglas parameters estimates by default in the function create.fleets.data
  #  when the stock is in biomass, therefore it is necessary to include 
  #  the historical gB data
  #       alpha=1, beta=1,q=Catches/(B+gB)*E
  
  stk2_gB.flq <- FLQuant(dimnames = list(age = 'all', year = first.yr:(proj.yr-1), unit = stk2.unit, 
                                         season = 1:ns, iter = 1:ni)) 
  
  
  first.yr <- 1990
  proj.yr <- 2009
  
  p <- 1
  r <- 1.185527
  K <- 89490
  stk2.biomass <- as.FLQuant(read.csv(file = 'data/stk2_biomass.csv'))
  for(ss in 1:ns){
    biomass <- as.vector(stk2.biomass[, as.character(first.yr:(proj.yr-1)),,ss,])
    bprod = biomass * (r/p) * (1 - (biomass/K)^p)
    stk2_gB.flq[,as.character(first.yr:(proj.yr-1)),,ss,] <- bprod
  }
  
  #         Projection
  #==============================================================================
  #         fleets: fl1
  #==============================================================================

    fl1_proj.avg.yrs           <- c(2003:2005)
    fl1.met1_proj.avg.yrs      <- c(2003:2005)   
    fl1.met2_proj.avg.yrs      <- c(2003:2005)   
    
    fl1.met1.stk1_proj.avg.yrs <- c(2003:2005)  
    fl1.met1.stk2_proj.avg.yrs <- c(2003:2005)
    
    fl1.met2.stk1_proj.avg.yrs <- c(2003:2005)  
    fl1.met2.stk2_proj.avg.yrs <- c(2003:2005)

  #==============================================================================
  #         fleets: fl2
  #==============================================================================
  
    fl2_proj.avg.yrs           <- c(2005,2007,2008)
    fl2.met1_proj.avg.yrs      <- c(2005,2007,2008)
    fl2.met2_proj.avg.yrs      <- c(2005,2007,2008)
    
    fl2.met1.stk1_proj.avg.yrs <- c(2005,2007,2008)
    fl2.met2.stk1_proj.avg.yrs <- c(2005,2007,2008)
    
    fl2.met1.stk2_proj.avg.yrs <- c(2005,2007,2008)
    fl2.met2.stk2_proj.avg.yrs <- c(2005,2007,2008)
    
    #  COBB douglas parameters
    
    #==============================================================================
    #              FLBEIA input object: fleets
    #==============================================================================
      
      fls.data <- list(fl1=ls(pattern="^fl1"),fl2=ls(pattern="^fl2")) 
      stks.data <- list(stk1=ls(pattern="^stk1"),stk2=ls(pattern="^stk2")) 
      
      fleets   <- create.fleets.data(yrs,ns,ni,fls.data,stks.data)
      for(fl in names(fleets)){
        
        fleets[[fl]]@fcost[] <- rnorm(prod(dim(fleets[[1]]@crewshare)), 5000,500)
        fleets[[fl]]@fcost[] <- rnorm(prod(dim(fleets[[2]]@crewshare)), 1000,10)

        fleets[[fl]]@crewshare[] <- rnorm(prod(dim(fleets[[1]]@crewshare)), .25,.025)
        fleets[[fl]]@metiers[[1]]@vcost[] <- rnorm(prod(dim(fleets[[1]]@crewshare)), 1000,100)
        fleets[[fl]]@metiers[[2]]@vcost[] <- rnorm(prod(dim(fleets[[2]]@crewshare)), 800,80)
        
        fleets[[fl]]@metiers[[1]]@catches[[1]]@price[] <- rnorm(length(c(fleets[[fl]][[1]][[1]]@price)),2500,250)
        fleets[[fl]]@metiers[[1]]@catches[[2]]@price[] <- rnorm(length(c(fleets[[fl]][[1]][[2]]@price)),2500,250)
        
        fleets[[fl]]@metiers[[2]]@catches[[1]]@price[] <- rnorm(prod(dim(fleets[[fl]][[2]][[1]]@price)),1800,180)
        fleets[[fl]]@metiers[[2]]@catches[[2]]@price[] <- rnorm(prod(dim(fleets[[fl]][[2]][[2]]@price)),800,80)
      
      }
      
      # # plot fleets
      # plotFLFleets(fleets,pdfnm='s0')   
      
      
#==============================================================================
#  Section 8:       SRs
#==============================================================================

  stk1_sr.model        <- 'bevholt'
  stk1_params.n        <- 2
  stk1_params.array    <- xtabs2(data~param+year+season+iter, 
                                 data=read.csv(file = 'data/stk1_params.csv'),exclude=NULL,na.action=na.pass)           
  stk1_params.name     <- c('a','b') 
  stk1_rec.flq         <- as.FLQuant(read.csv(file = 'data/stk1_rec.csv'))             
  stk1_ssb.flq         <- as.FLQuant(read.csv(file = 'data/stk1_ssb.csv'))             
  stk1_uncertainty.flq <- as.FLQuant(read.csv(file = 'data/stk1_uncertainty.csv'))     
  stk1_proportion.flq  <- as.FLQuant(read.csv(file = 'data/stk1_proportion.csv'))      
  stk1_prop.avg.yrs    <- ac(2006:2008)
  stk1_timelag.matrix  <- matrix(c(0,1),2,ns, dimnames = list(c('year', 'season'),season = 1:ns)) 
  
  #==============================================================================
  #              FLBEIA input object: SRs
  #==============================================================================
  
    stks.data <- list(stk1=ls(pattern="^stk1"))
    
    SRs <- create.SRs.data(yrs,ns,ni,stks.data)
    
    
#==============================================================================
#  Section 9:       BDs
#==============================================================================

  stk2_bd.model        <- 'PellaTom'
  stk2_params.name     <- c('K','p','r') 
  stk2_params.array    <- xtabs2(data~param+year+season+iter, 
                                 data=read.csv(file = 'data/stk2_params.csv'),exclude=NULL,na.action=na.pass)              
  stk2_biomass.flq     <- as.FLQuant(read.csv(file = 'data/stk2_biomass.csv'))            
  stk2_catch.flq       <- as.FLQuant(read.csv(file = 'data/stk2_catch.csv'))               
  stk2_uncertainty.flq <- as.FLQuant(read.csv(file = 'data/stk2_uncertainty.csv'))          
  stk2_alpha <- array(1,dim=c(length(first.yr:last.yr),ns,ni))
  #    stk2_gB.flq     already defined in line 327
  # Calculate gB in the biomass dynamic population
  
  #==============================================================================
  #              FLBEIA input object: BDs
  #==============================================================================
  
    stks.data <- list(stk2=ls(pattern="^stk2")) 
    BDs       <- create.BDs.data(yrs,ns,ni,stks.data)
    
    BDs$stk2@alpha <- array((BDs$stk2@params['p',,,]/BDs$stk2@params['r',,,]+1)^(1/ BDs$stk2@params['p',,,]), dim = c(length(first.yr:last.yr),4,1))
    
    for(y in 2:19){
      for(s in 1:4){
        BDs[[1]] <- FLBEIA:::BDsim(BDs[[1]], year = y, season = s)
      }
    }
     
     
#==============================================================================
#  Section 10:       advice:TAC/TAE/quota.share
#==============================================================================

  stk1_advice.TAC.flq <- as.FLQuant(read.csv(file = 'data/stk1.tac.csv'))
  stk1_advice.avg.yrs <- c(2006:2008)
  
  stk2_advice.TAC.flq <- as.FLQuant(read.csv(file = 'data/stk2.tac.csv'))
  stk2_advice.avg.yrs <- c(2006:2008)
  
  #==============================================================================
  #              FLBEIA input object: advice
  #==============================================================================
  
    stks.data <- list(stk1=ls(pattern="^stk1"),stk2=ls(pattern="^stk2"))
    
    advice <- create.advice.data(yrs,ns,ni,stks.data,fleets)
    
    
#==============================================================================
#  Section 11:       main.ctrl
#==============================================================================

  main.ctrl           <- list()
  main.ctrl$sim.years <- c(initial = proj.yr, final = last.yr)
  
  
#==============================================================================
#  Section 12:       biols.ctrl
#==============================================================================

  growth.model <- c('ASPG','BDPG')
  biols.ctrl <- create.biols.ctrl (stksnames=stks,growth.model=c('ASPG','BDPG'))
  

#==============================================================================
#  Section 13:       fleets.ctrl
#==============================================================================

  n.fls.stks       <- c(2,2)
  fls.stksnames    <- rep(c('stk1','stk2'),2)
  effort.models    <- c('fixedEffort','SMFB')
  effort.restr.fl2 <- 'min'
  restriction.fl2  <- 'catch'
  catch.models     <- rep(c("CobbDouglasAge", "CobbDouglasBio"),2)
  capital.models   <- rep('fixedCapital',2)
  price.models     <- NULL                            
  
  flq.stk1 <- FLQuant(dimnames = list(age = 'all', year = first.yr:last.yr, unit = stk1.unit, 
                                      season = 1:ns, iter = 1:ni)) 
  fleets.ctrl <- create.fleets.ctrl(fls=fls,n.fls.stks=n.fls.stks,fls.stksnames=fls.stksnames,
                                    effort.models= effort.models, catch.models=catch.models,
                                    capital.models=capital.models, price.models=price.models,flq=flq.stk1,
                                    effort.restr.fl2 = effort.restr.fl2, restriction.fl2 = restriction.fl2)
  
  fleets.ctrl$fl1$restriction <- "landings"
  fleets.ctrl$fl2$restriction <- "landings"
  
    
#==============================================================================
#  Section 14:       advice.ctrl
#==============================================================================

  HCR.models <- c('IcesHCR','annualTAC') 
  Ftarget.stk2 <- 0.3
  ref.pts.stk1 <- matrix(c(0.24, 50000, 150000), 3,ni, dimnames = list(c('Fmsy', 'Blim', 'Btrigger'), 1:ni))
  
  advice.ctrl <- create.advice.ctrl(stksnames = stks, HCR.models =  HCR.models,
                                    ref.pts.stk1 = ref.pts.stk1, Ftarget.stk2 = Ftarget.stk2,first.yr=first.yr,last.yr=last.yr)
  
  advice.ctrl[['stk1']][['sr']]            <- list()
  advice.ctrl[['stk1']][['sr']][['model']] <- 'segreg'
  advice.ctrl[['stk1']][['sr']][['years']] <- c(y.rm = 2, num.years = 10)
  advice.ctrl$stk1$AdvCatch <- rep(TRUE,length(first.yr:last.yr))   #TRUE advice in catches, FALSE advice in landings
  names(advice.ctrl$stk1$AdvCatch) <- as.character((first.yr:last.yr))
  

#==============================================================================
#  Section 15:       assess.ctrl
#==============================================================================

  assess.models <- rep('NoAssessment',2)
  
  assess.ctrl <- create.assess.ctrl(stksnames = stks, assess.models = assess.models)
  

#==============================================================================
#  Section 16:       obs.ctrl
#==============================================================================

  stkObs.models <- rep('perfectObs',2)
  
  obs.ctrl <- create.obs.ctrl(stksnames = stks,  stkObs.models = stkObs.models)
  

#==============================================================================
#  Section 17:       covars
#==============================================================================

  cv_mean_value <- c( FuelCost = 46, CapitalCost = 4519.06, Salaries = 0, InvestShare = 0.2, NumbVessels = 228.33, 
                      MaxDays = 228, w1 = 0.03, w2 = 0.03, EmploymentPerVessel = 2)
  
  covars <- vector("list",9)
  names(covars) <- c("FuelCost","CapitalCost","Salaries", "InvestShare","NumbVessels","MaxDays",
                     "w1","w2","EmploymentPerVessel")
  flq <- FLQuant(rnorm(length(fls)*length(first.yr:last.yr)*ns*ni, 1000,100), 
                 dimnames = list(fleets = fls, year = first.yr:last.yr, unit = stk1.unit, season = 1:ns, iter = 1:ni)) 
  
  for(cv in names(covars)){ 
    covars[[cv]] <- flq
    dimnames(covars[[cv]])[[1]] <- c('fl1', 'fl2') 
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
    covars.ctrl[[cv]][['process.model']] <- 'fixedCovar'
  }
  
  
#==============================================================================
#  Section 18:       FLBEIA input objects
#==============================================================================

  indices <- NULL
  

#==============================================================================
#       Rename
#==============================================================================

  multiAdv   <- advice      
  multiAdvC  <- advice.ctrl 
  multiAssC  <- assess.ctrl
  multiBD    <- BDs
  multiBio   <-  biols       
  multiBioC  <- biols.ctrl  
  multiCv    <- covars
  multiCvC   <- covars.ctrl
  multiFl    <- fleets      
  multiFlC   <- fleets.ctrl 
  multiMainC <- main.ctrl  
  multiObsC  <-  obs.ctrl
  multiSR    <- SRs
  
  
#==============================================================================
#        FLBEIA input objects
#==============================================================================
  
  save( multiAdv, multiAdvC, multiAssC, multiBD, multiBio, multiBioC, multiCv, multiCvC, 
        multiFl, multiFlC, multiMainC, multiObsC, multiSR, 
        file = '../../data/multi.RData')
    
  
# #==============================================================================
# #  Section 19:       Run FLBEIA
# #==============================================================================
# 
#   s2 <- FLBEIA(biols = multiBio,   # FLBiols object with 2 FLBiol element for stk1.
#                SRs = multiSR,      # A list with 1 FLSRSim object for stk1.
#                BDs = multiBD,      # A list with 1 FLBDSim object for stk2.
#                fleets = multiFl,   # FLFleets object with on fleet.
#                covars = multiCv,   # covars not used
#                indices = NULL,     # indices not used 
#                advice = multiAdv,  # A list with two elements 'TAC' and 'quota.share'
#                main.ctrl = multiMainC,  # A list with one element to define the start and end of the simulation.
#                biols.ctrl = multiBioC,  # A list with one element to select the model to simulate the stock dynamics.
#                fleets.ctrl = multiFlC,  # A list with several elements to select fleet dynamic models and store additional parameters.
#                covars.ctrl = multiCvC,  # covars control not used 
#                obs.ctrl = multiObsC,    # A list with one element to define how the stock observed ("PerfectObs").
#                assess.ctrl = multiAssC, # A list with one element to define how the stock assessment model used ("NoAssessment").
#                advice.ctrl = multiAdvC) # A list with one element to define how the TAC advice is obtained ("IcesHCR").
#   
#   multiRes <- s2
#   save(multiRes,file='../../data/multiRes.RData')
#   
#   
# #==============================================================================
# #  Section 20:       Plot
# #==============================================================================
# 
#   # Names of the object returned by FLBEIA
#   names(s2)
#   
#   # The default plot for FLBiol defined in FLCore
#   plot(s2$biols[[1]])
#   
#   # Create summary data frames (biological, economic, and catch)
#   
#   s2_sum      <- bioSum(s2)
#   s2_flt      <- fltSum(s2)
#   
#   s2_flt      <- fltSum(s2, byyear = FALSE)
#   
#   s2_fltStk   <- fltStkSum(s2)
#   
#   
#   # Create several plots and save them in the working directory using 'pdf' format and 
#   # 's2' suffix in the name.
#   
#   plotFLBiols(s2$biols, pdfnm='s2')
#   plotFLFleets(s2$fleets,pdfnm='s2')
#   plotfltStkSum(s2,pdfnm='s2') 
#   plotEco(s2$fleets,pdfnm='s2')
#   
#   ## End(Not run)
    
  