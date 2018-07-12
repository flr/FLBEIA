#-------------------------------------------------------------------------------  
# create.biol.arrays function: 
# function to generate an FLFleets object given inputs as arrays
# 
# Created: Dorleta Garcia -  2018-04-11
# Changed: 2018-07-04 10:58:14 (ssanchez)
#------------------------------------------------------------------------------- 

# Conditioning_create.biol.arrays.r - function to calculate indicesB and indicesP (given some information on growth, periods, catch and an FLPar object)
# FLBEIA/R/Conditioning_create.biol.arrays.r

# Copyright: AZTI, 2018
# Author: Dorleta Garcia & Sonia Sanchez (AZTI) (<dgarcia@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# create.biol.arrays {{{

#' @title Function to generate an FLBiol object given inputs as arrays
#'
#' @description This function generates an FLBiol object, given the data inputs as arrays. 
#'              Supported formats are Excel (xls and xlsx) and R format (RData).
#'
#' @name create.biol.arrays
#' @rdname create.biol.arrays
#' @aliases create.biol.arrays
#'
#' @param filename A character vector with the name of the files containing the stock data.
#'                 Supported formats are Excel (xls and xlsx) and R format (RData).
#'                 In case of using R format, the information must be stored in \code{data} object (consisting in a list with the different elements). 
#'                 The following information is compulsory: abundances in numbers at age (n), mean weight at age (wt), maturity (mat), 
#'                 natural mortality (m), moment of the year when spawning occurs in percentage (spwn), fishing mortality at age (f) and 
#'                 catch in numbers at age (caa).
#'                 For the rest of information, if not provided, default values are set. For example, fecundity (fec) is set to 1,
#'                 landings and discard in numbers at age (laa and daa) are set to cca and 0, respectively. 
#'                 Finally for weights, if missing, weights at age for landings (wl) and discards (wd) are set to the weights in the population and 
#'                 weights at age for catch (wc) are set to the weighted mean of the weights of landings and discards.
#' @param name     A character (optional) with the name of the stock.
#' @param ages     A numeric vector with the age classes of stock.
#' @param fbar     A numeric vector with the age range (min,max) to be used for estimating average fishing mortality.
#' @param hist.yrs A vector with the historical years.
#' @param sim.yrs  A vector with the simulation years.
#' @param mean.yrs A vector with the years used to compute the mean to condition the parameters in the projection period.
#' @param excel    Logical. \code{TRUE} (default), if the data is provided in an Excel file and 
#'                 \code{FALSE}, if an RData object is used instead.
#' @param unit     A list with the units of the different elements included in \code{filename}. Unitless objects must be set to '' or character(1).
#'                 This parameter is only required if \code{excel==FALSE}. When using Excell files the units are taken from the first row and column (cell A1) of each sheet. 
#'                 If the cell is empty then units are set to NA, in case of an unitless object then 1 must be inputed into cell A1.
#' 
#' @return An \code{FLBiol}. 
#'
#' @author Dorleta Garcia & Sonia Sanchez.
#' @seealso \link{FLBiol}, \link{create.fleets.arrays}
#' @keywords create.biol.arrays
#'
#'  
# @examples 
# 
# # still missing an example
# 



create.biol.arrays <- function(filename, name = NA, ages, hist.yrs , sim.yrs, fbar = NULL, mean.yrs, excel = TRUE, unit = list()){
  
  ages     <- ac(ages)
  hist.yrs <- ac(hist.yrs)
  sim.yrs  <- ac(sim.yrs)
  mean.yrs  <- ac(mean.yrs)
  
  if (sum(!mean.yrs %in% hist.yrs)>0) stop('mean.yrs must be taken from hist.yrs')
  
  sheets <- c('n', 'wt', 'mat', 'fec', 'm', 'spwn')
  
  data <- vector('list', 6)
  names(data) <- sheets
  
  if (length(ages)==1) {
    if (ages!='all') {
      ages <- 'all'
      warning("ages has been renamed to 'all'")
    }
    if (!is.null(fbar)) {
      fbar <- NULL
      warning("fbar has been set to c(1,1)")
    }
  } else if (is.null(fbar)) {
    fbar <- c(ages[1], ages[length(ages)])
  } else {
    # check format
    if (length(fbar)>2)
      stop('fbar must be a vector of the form: c(minfbar,maxfbar)')
    # check ranges
    if ( fbar[1]<as.numeric(ages[1])) stop(paste('minfbar must be >=',ages[1]))
    if ( fbar[2]>as.numeric(ages[length(ages)])) stop(paste('maxfbar must be <=',ages[length(ages)]))
  }
  
  yrs <- hist.yrs[1]:sim.yrs[length(sim.yrs)]
  
  nage  <- length(ages)
  nyear <- length(yrs)
  
  if(excel == TRUE){
    if (!is.null(names(unit))) 
      stop('units must be set in the Excel file')
    wb <- loadWorkbook(filename, create = FALSE)
    wb_sheets <- getSheets(wb)
    # check that all required sheets are available
    if ( any(!sheets[sheets!="fec"] %in% wb_sheets))
      stop(paste("Sheets: ", paste(sheets[sheets!="fec"], collapse = ", "), " are required in file: '", filename, "'", sep=''))
    for(sl in sheets)  {
      if (sl=='fec' & (!sl %in% wb_sheets)) { # if missing fec --> set equal to 1
        data[[sl]] <- data[['mat']]*0+1
        unit[[sl]] <- ''
        next()
      }
      # check ages
      aa <- readWorksheet(wb, sheet = sl, header = FALSE, startRow = 2, startCol = 1, endCol = 1)$Col1
      if (length(ages)!=length(aa)) {
        stop(paste("check age range in sheet '",sl,"' as it is different from 'ages'"),sep='')
      } else if (length(ages)>1 & sum(ages!=aa)>0) 
        stop(paste("ages in sheet '",sl,"' are different from ",ages[1],":",ages[length(ages)],sep=''))
      # check years
      yy <- unlist(readWorksheet(wb, sheet = sl, header = FALSE, startRow = 1, startCol = 2, 
                                                          endRow = 1, endCol = nyear + 1))
      if (sum(hist.yrs!=yy)>0) 
        stop(paste("years in sheet '",sl,"' are different from ",hist.yrs[1],":",hist.yrs[length(hist.yrs)],sep=''))
      data[[sl]] <- as.matrix(readWorksheet(wb, sheet = sl, header = TRUE, startRow = 1, startCol = 2, 
                                                                            endRow = nage + 1, endCol = nyear + 1))
      unit[[sl]] <- readWorksheet(wb, sheet = sl, header = FALSE, startRow = 1, startCol = 1, endRow = 1, endCol = 1)$Col1
    }
    nit  <- 1
    unit <- lapply( unit, function(x) ifelse( is.na(x), 'NA', ifelse( x==1 | x=='1', '', as.character(x))))
  }
  else{
    data <- loadToEnv(filename)[["data"]]
    nit <- ifelse(is.na(dim(data$n)[3]), 1, dim(data$n)[3])
    if (is.null(names(unit)))
      warning('Please remember to set the units for the different slots!')
  }
  
  
  flq <- FLQuant(dim = c(length(ages), length(yrs), 1,1,1,nit), dimnames = list(age = ages, year = yrs, iter = 1:nit))
  
  if (length(ages)==1) {
    res <- FLBiol(name = name, 
                  desc = paste('data imported from', filename), 
                  range = c(min = NA, max = NA, plusgroup = NA,  
                            minyear = as.numeric(hist.yrs[1]), maxyear = as.numeric(sim.yrs[length(sim.yrs)]), 
                            minfbar = 1, maxfbar = 1),
                  spwn = flq)
  } else {
    res <- FLBiol(name = name, 
                  desc = paste('data imported from', filename), 
                  range = c(min = as.numeric(ages[1]), max = as.numeric(ages[length(ages)]), plusgroup = as.numeric(ages[length(ages)]),  
                            minyear = as.numeric(hist.yrs[1]), maxyear = as.numeric(sim.yrs[length(sim.yrs)]), 
                            minfbar = fbar[1], maxfbar = fbar[2]),
                  spwn = flq)
  }
  
  res@n[,hist.yrs]  <- data$n
  res@m[,hist.yrs]  <- data$m
  res@wt[,hist.yrs] <- data$wt
  res@spwn[,hist.yrs]  <- data$spwn
  for(sl in c('n', 'wt', 'm', 'spwn')) {
    units(res)[[sl]] <- unit[[sl]]
  }
  res@rec <- predictModel(n = res@n, model = ~ n[1,])
  mat(res)[,hist.yrs] <- data$mat #res@mat$mat[,hist.yrs] <- data$mat
  if (!is.null(unit[['mat']]))
    units(mat(res)) <- unit[['mat']]
  fec(res)[,hist.yrs] <- data$fec #res@fec$fec[,hist.yrs] <- data$fec
  if (!is.null(unit[['fec']])) 
    units(fec(res)) <- unit[['fec']]

  
  # projection
  res@m[,sim.yrs]       <- yearMeans(res@m[,mean.yrs])
  res@wt[,sim.yrs]      <- yearMeans(res@wt[,mean.yrs])
  res@spwn[,sim.yrs]    <- yearMeans(res@spwn[,mean.yrs])
  mat(res)[,sim.yrs] <- yearMeans(res@mat$mat[,mean.yrs]) #res@mat$mat[,sim.yrs] <- yearMeans(res@mat$mat[,mean.yrs])
  fec(res)[,sim.yrs] <- yearMeans(res@fec$fec[,mean.yrs]) #res@fec$fec[,sim.yrs] <- yearMeans(res@fec$fec[,mean.yrs])
  
  return(res)
         
}

