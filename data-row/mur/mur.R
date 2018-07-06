#-------------------------------------------------------------------------------  
# mur dataset: 
# code to generate data in mur.RData
# 
# Created: Agurtzane Urtizberea -  2016-02-09
# Changed: 2018-07-06 11:15:22 (ssanchez)
#------------------------------------------------------------------------------- 

# mur.r - code to generate data in mur.RData
# FLBEIA/data-row/mur/mur.r

# Copyright: AZTI, 2018
# Author: Dorleta Garcia and Sonia Sanchez (AZTI) (<flbeia@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


#==============================================================================
# Load the data
#==============================================================================

  catch <- read.csv("mur_catch.txt", sep = " ") 
  names(catch) <- c("year", paste('area', c("9a",  "5b", "6ab",  "celtic",  "bob", "channel", "total"), sep = "_"))             
  
  catch_long <- reshape(catch, direction = 'long', varying = names(catch)[2:8], sep = "_")
  names(catch_long)[3:4] <- c('area', 'catch') 
  
  evhoe <- read.csv("Biomass_IndexEvohe_BoB_CS.csv", sep = ";")[,5:8]
  names(evhoe) <- c('year', 'biomass', 'std','cv')
  
  # In 1999 France did not report any data!!
  
  catch[,'year'] <- as.numeric(as.character(catch[,'year']))
  catch[25, 'year'] <- 1999
  catch[25, 'area_total'] <- NA


#==============================================================================
# Save the data
#==============================================================================

  save( catch, evhoe, file='../../data/mur.RData')

