#-------------------------------------------------------------------------------  
# create.ecoData function: 
# function to generate an covars object with economic indicators 
# 
# Created: Dorleta Garcia -  2018-11-22
#------------------------------------------------------------------------------- 

# Copyright: AZTI, 2018
# Author: Dorleta Garcia (<dgarcia@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# create.ecoData {{{

#' @title Function to generate the covar object neccesary to compute economic indicators and fills in economic slots in FLFleets object 
#' 
#'
#' @description This function generates an covar object with the economic indicators and fills in vcost, fcost, capacity and crewshare 
#'             slots in FLFleets object. The data is given in excel File (xls and xlsx).
#'
#' @name create.ecoData
#' @rdname create.ecoData
#' @aliases create.ecoData
#'
#' @param file            An excel file with the economic data. The names used for the sheets and the columns must be the same used to name fleets and metiers in fltObj. However, the fleets in the excel file can be a subset of the fleets in the fleetsObj, i.e., is not neccesary to provide economic data for all the fleets. The order of metiers in the columns data must correspond with the name of the metiers used in the FLFleet object.
#' @param fltObj         An FLFleets object  with the structure of the fleet and which may contain historical data. 
#' @param hist.yrs        A vector with the historical years.
#' @param sim.yrs         A vector with the simulation years.
#' @param mean.yrs        A vector with the years used to compute the mean to condition the parameters in the projection period.
#' 
#' @return An /code{FLFleetsExt}.
#'
#' 
#' @author Dorleta Garcia & Sonia Sanchez.
#' @seealso \code{\link{FLFleetsExt}}, \code{\link{create.biol.arrays}}
#' @keywords create.fleets.arrays
#'
#'  
# @examples 
# 
# # still missing an example
# 


# 1. fuel cost
# 2. crew share
# 3. other variable cost
# 4. fixed costs
# 5. capital value
# 6. fixed salarie
# 7. maximum effort
# 8. crew
# 9. depreciation
# 10. vessels
# 11. new vessel
# 12. investment share
# 13. w1
# 14. w2


create.ecoData <- function(file, fltObj, hist.yrs, mean.yrs, sim.yrs){
  
  flObj <- unlist(fltObj)
  
  wb <- loadWorkbook(file, create = FALSE)
  sheets <- getSheets(wb)
  sheets <- sheets[!(sheets == 'readme')]

  fltnms <- names(fltObj)
  if(!all(sheets %in% fltnms)) stop('Some of the fleets in the excel file do not correspond with the fleets in the FLFleetsObj (fltObj)')
  
  # Structure of the covars Object
  flq <- FLQuant(dimnames = list(fleet = sheets, year = c(hist.yrs, sim.yrs)))
  
  covars <- lapply(list("FuelCost", "CapitalValue", "Salaries", "InvestShare", "NumbVessels", "MaxDays",
                     "w1","w2", "EmploymentPerVessel", "NewVesselCost"), function(x) return(flq))
  names(covars) <- c("FuelCost", "CapitalValue", "Salaries", "InvestShare", "NumbVessels", "MaxDays",
                     "w1","w2", "EmploymentPerVessel", "NewVesselCost")
    
  # Fill the flObj
  for(fl in sheets){
    dat <- readWorksheet(wb, sheet = fl, header = TRUE, startRow = 1, startCol = 1, endRow =  15)
    dat$indicator <- tolower(dat$indicator) 
    names(dat)[-(1:3)] <- names(fltObj[[fl]]@metiers)
      
    # check if all the metiers are available, if not stop the function.
    if(!(identical(names(fltObj[[fl]]@metiers), names(dat)[-(1:3)]))) stop(cat('The names of the metiers in fleet ', fl, ' do not corresponds with the names used in the excel file.\n'))
    
    # We order the data by year and indicator, in the code below it is assumed that there is the same number of years per indicator.
    dat <- dat[order(dat$indicator),]
    dat <- dat[order(dat$year),]
    for(k in 3:dim(dat)[2]) dat[,k] <- as.numeric(dat[,k])
    
    # fleet level
    fcost     <- subset(dat, indicator == 'fixed costs')
    crewshare <- subset(dat, indicator == 'crew share')
    nvessels  <- subset(dat, indicator == 'vessels')
    maxef     <- subset(dat, indicator == 'maximum effort')
    
    fltObj[[fl]]@fcost[,ac(fcost$year)] <- c(fcost[,'fleet'])
    fltObj[[fl]]@crewshare[,ac(crewshare$year)] <- c(crewshare[,'fleet'])
    fltObj[[fl]]@capacity[,ac(nvessels$year)] <- c(nvessels[,'fleet']*maxef[,'fleet'])
    
    # metier level
    vcost     <- subset(dat, indicator == 'fuel cost')
    otvcost   <- subset(dat, indicator == 'other variable cost')
    
    for(mt in names(fltObj[[fl]]@metiers))  fltObj[[fl]]@metiers[[mt]]@vcost[,ac(vcost$year)] <- vcost[,mt] + otvcost[,mt]
    
    # Fill in the covars object for fl fleet.  
    fuC <- subset(dat, indicator == 'fuel cost')
    caV <- subset(dat, indicator == 'capital value')
    fiS <- subset(dat, indicator == 'fixed salarie')
    maE <- subset(dat, indicator == 'maximum effort')
    dep <- subset(dat, indicator == 'depreciation')
    ves <- subset(dat, indicator == 'vessels')
    nev <- subset(dat, indicator == 'new vessel')
    inv <- subset(dat, indicator == 'investshare')
    emp <- subset(dat, indicator == 'employees')
    w1  <- subset(dat, indicator == 'w1')
    w2  <- subset(dat, indicator == 'w2')
  
    # FuelCost is a weighted mean of the fuelcost in each metier.
    covars[["FuelCost"]][fl,ac(fuC$year)] <- apply(matrix(sapply(fltObj[[fl]]@metiers, 
                                                          function(x) c(x@effshare[,ac(fuC$year)])*c(fuC[,x@name])),length(fuC$year),length(fltObj[[fl]]@metiers)),
                                                   1,sum)
    covars[["CapitalValue"]][fl,ac(caV$year)]        <- caV[,'fleet']
    covars[["Salaries"]][fl,ac(fiS$year)]            <- fiS[,'fleet']
    covars[["InvestShare"]][fl,ac(inv$year)]         <- inv[,'fleet']
    covars[["NumbVessels"]][fl,ac(ves$year)]         <- ves[,'fleet']
    covars[["MaxDays"]][fl,ac(maE$year)]             <- maE[,'fleet']
    covars[["w1"]][fl,ac(w1$year)]                   <- w1[,'fleet']
    covars[["w2"]][fl,ac(w2$year)]                   <- w2[,'fleet']
    covars[["NewVesselCost"]][fl,ac(nev$year)]       <- nev[,'fleet']
    covars[["EmploymentPerVessel"]][fl,ac(emp$year)] <- emp[,'fleet']
  }
    
    #---------------------------------------------------------------------------
    # Values in the projection
    #---------------------------------------------------------------------------
    # FLFleetsExt object
    for(fl in sheets){
      for(sl in c('crewshare', 'fcost', 'capacity')){
          slot(fltObj[[fl]],sl)[, ac(sim.yrs)] <-  yearMeans(slot(fltObj[[fl]],sl)[, ac(mean.yrs)])
  
          if(sl == 'capacity' & any(is.na(c( slot(fltObj[[fl]],sl)[, ac(sim.yrs)] )))){
            fltObj[[fl]]@capacity[, is.na(fltObj[[fl]]@capacity[, ac(sim.yrs)])] <- 1e12
            cat('Warning: Capacity slot in fleet ', fl, ' has been set to 1e12 to avoid problems in the projection.\n')
          }
      }
      for(mt in names(fltObj[[fl]]@metiers)){
         slot(fltObj[[fl]]@metiers[[mt]], 'vcost')[, ac(sim.yrs)] <-  yearMeans(slot(fltObj[[fl]]@metiers[[mt]],'vcost')[, ac(mean.yrs)])
          
      }
    }
    #covars object
    for(cv in c("FuelCost", "CapitalValue", "Salaries", "InvestShare", "NumbVessels",         
                "MaxDays", "w1", "w2", "EmploymentPerVessel")){
      covars[[cv]][, ac(sim.yrs)] <-  yearMeans(covars[[cv]][, ac(mean.yrs)])
    }
  

  


  return(list(fleets = FLFleetsExt(fltObj), covars = covars))
}


#
# TEST
#
# load("~/OneDrive - AZTI/BoB/02_MixedFisheries/data/FLR_Objs.RData")
# library(XLConnectJars)
# library(XLConnect)
# library(FLBEIA)
# res <- create.ecoData('C:/use/OneDrive - AZTI/BoB/03_Economics/input/eco/economic_data.xlsx', fleets, 2005:2017, 2017, 2018:2025)




