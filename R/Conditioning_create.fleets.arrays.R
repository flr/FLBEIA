#-------------------------------------------------------------------------------  
# create.fleets.arrays function: 
# function to generate an FLFleets object given inputs as arrays
# 
# Created: Dorleta Garcia -  2018-04-11
# Changed: 2018-07-04 10:58:14 (ssanchez)
#------------------------------------------------------------------------------- 

# Conditioning_create.fleets.arrays.r - function to calculate indicesB and indicesP (given some information on growth, periods, catch and an FLPar object)
# FLBEIA/R/Conditioning_create.fleets.arrays.r

# Copyright: AZTI, 2018
# Author: Dorleta Garcia & Sonia Sanchez (AZTI) (<dgarcia@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# create.fleets.arrays {{{

#' @title Function to generate an FLFleetsExt object given inputs as arrays
#'
#' @description This function generates an FLFleetsExt object, given the data inputs as arrays. 
#'              Supported formats are Excel (xls and xlsx) and R format (RData).
#'
#' @name create.fleets.arrays
#' @rdname create.fleets.arrays
#' @aliases create.fleets.arrays
#'
#' @param stk_objs        A character vector with the names of the files containing the stocks data. See \link{create.biols.arrays} for more detail.
#'                        Supported format is only Excel (xls and xlsx), each stock can be in different format.
#' @param caa_objs        A character vector with the names of the files containing the catch at age data (in numbers), both for landings and discards.
#'                        Supported formats are Excel (xls and xlsx) and R format (RData), each file can be in different format.
#'                        The number of required files depend of the value of \code{caaOpt} argument:
#'                        \itemize{ 
#'                           \item{\code{caaOpt} = 1 or 2   }{=> one per stock and fleet. Named vector stknm_fltnm}
#'                           \item{\code{caaOpt} = 3, 4 or 5}{=> one per stock. Named vector stknm}
#'                        }
#'                        If NULL, the function looks for "caa_stknm_flnm.xlsx" in \code{caaOpt} = 1 or 2 and "caa_stknm.xlsx" in \code{caaOpt} = 3 ,4 or 5.
#' @param caa_objs_path   A character vector with the \code{caa_objs} file path.
#' @param update_price    Logical. If \code{TRUE} (default), prices must be provided and \code{price_objs} and \code{price_objs_path} arguments are required.
#' @param price_objs      A character vector with the names of the files containing the price at age data.
#'                        Supported formats are Excel (xls and xlsx) and R format (RData), each prices file can be in different format.
#'                        The number of required files depend of the value of \code{caaOpt} argument:
#'                        \itemize{ 
#'                           \item{\code{caaOpt} = 1 or 2   }{=> one per stock and fleet. Named vector stknm_fltnm.}
#'                           \item{\code{caaOpt} = 3, 4 or 5}{=> one per stock. Named vector stknm.}
#'                        }
#' @param price_objs_path A character vector with the price_objs file path.
#' @param catch_obj       A character vector with the names of the files containing the catch data in Fcube format.
#'                        Supported formats are Excel (xls and xlsx) and R format (RData) and required columns are 'year', 'fleet', 'metier', 'stock', 'category' and 'catch'.
#' @param effort_obj      A character vector with the names of the files containing the effort data in Fcube format.
#'                        Supported formats are Excel (xls and xlsx) and R format (RData) and required columns are 'year', 'fleet', 'metier' and 'effort'.
#'                        Both \code{catch_obj} and \code{effort_obj} must be in the same format.
# @param eco_obj         A character vector with the names of the files containing the economic data in Fcube format.
#                        Supported formats are Excel (xls and xlsx) and R format (RData) and required columns are 'year', 'fleet', 'metier' and 'effort'.
#                        Both \code{catch_obj} and \code{effort_obj} must be in the same format.
# @param stk_obj         An FLBiols object (optional), with one FLBiol per stock.
#                        If this object is provided, then +++
#' @param flt_obj         An FLFleets object (optional) with the structure of the fleet and which may contain historical data. 
#'                        If this object is provided, then the arguments \code{stk_nms}, \code{flt_nms}, \code{flt_mt_nms} and \code{ages} will not be used
#'                        and \code{excel} argument must be set to \code{FALSE}.
#' @param stk_nms         A character vector (optional) with the name of all the stocks caugth by the different fleets.
#' @param flt_nms         A character vector with the name of the fleets.
#' @param flt_mt_nms      A list with one element per fleet. In turn, each element is a character vector with the names of the metiers in the corresponding fleet.
#' @param flt_mt_stk_nms  A list with one element per fleet and metier. In turn, each element is a character vector with the names of the stocks in the corresponding fleet and metier.
#' @param ages            A list with one element per stock, with the age classes of the stock.
#' @param hist.yrs        A vector with the historical years.
#' @param sim.yrs         A vector with the simulation years.
#' @param mean.yrs        A vector with the years used to compute the mean to condition the parameters in the projection period.
#' @param new_hist.yrs    A vector with the years from input files that will be used to condition the parameters in the historic years.
#'                        If a value is not provided, the it is set equal to \code{hist.yrs}.
#' @param caa_flt_mt_correspondences An Excel file name. This file must contain one sheet per stock, with the correspondences between the fleet segments used in \code{caaa_obj} data 
#'                                   and the fleet metier segmentation used in the analysis. If the file does not exist, it is supposed that the caa data is given by fleet and metier. 
#' @param paa_flt_mt_correspondences An Excel file name. This file must contain information on prices correspondences, with same format and requirements as \code{caa_flt_mt_correspondences} argument.
#' @param caaOpt          A code number to determine the way in wich catch at age data are provided.
#'                        The option to be used depends on the data availabiltiy, from data rich to data-poor and the following codes are available:
#'                        \itemize{ 
#'                           \item{1}{If catch at age data is available at métier level for all the métiers.}
#'                           \item{2}{If catch at age data is only available at fleet level.}
#'                           \item{3}{If catch at age data is disaggregated but the segments do not correspond exactly with the métiers/fleets considered in the case study.}
#'                           \item{4}{If catch at age data is only available at stock level.}
#'                           \item{5}{If we want to use the data available previously in the \code{FLCatch} objects from \code{flt_obj} to derive catch profiles at age 
#'                                    and then apply \code{caaOpt==3} using only one fleet segment, fseg, which represents all the fleets and métiers.
#'                                    Note: This approach could lead to a different total catch at age profile derived from the fleets to those in the stocks.}
#'                        }
#' @param priceOpt        A code number to determine the way in wich price at age data are provided.
#'                        The option to be used depends on the data availabiltiy, from data rich to data-poor and the following codes are available:
#'                        \itemize{ 
#'                           \item{1}{If price data is available at métier level for all the métiers.}
#'                           \item{2}{If price data is only available at fleet level.}
#'                           \item{3}{If price data is disaggregated but the segments do not correspond exactly with the métiers/fleets considered in the case study.}
#'                           \item{4}{If price data is only available at stock level.}
#'                        }
#' @param excel           Logical. If \code{TRUE} (default), the data in the Excel file is used to create the stucture of the FLFleets object; whereas 
#'                        if \code{FALSE}, the \code{flt_obj} object is used instead.
#' 
#' @return An \code{FLFleetsExt}.
#'
#' 
#' @author Dorleta Garcia & Sonia Sanchez.
#' @seealso \link{FLFleetsExt}, \link{create.biols.arrays}
#' @keywords create.fleets.arrays
#'
#'  
# @examples 
# 
# # still missing an example
# 

create.fleets.arrays <- function(stk_objs, caa_objs, caa_objs_path, update_price = TRUE, price_objs, price_objs_path, 
                                 catch_obj, effort_obj, flt_obj = NULL, # stk_obj= NULL, eco_obj=NULL,
                                 stk_nms = NA, flt_nms, flt_mt_nms, flt_mt_stk_nms, 
                                 ages, hist.yrs, sim.yrs, mean.yrs, new_hist.yrs = hist.yrs,
                                 caa_flt_mt_correspondences = NULL, paa_flt_mt_correspondences = NULL, caaOpt, priceOpt, excel = TRUE){
  
  ages_stk <- lapply(ages, function(x) ac(x))
  hist.yrs <- ac(hist.yrs)
  sim.yrs  <- ac(sim.yrs)
  mean.yrs <- ac(mean.yrs)
  
  stks <- names(ages_stk)
  
  # stk_data <- caa_data <- list()
  
  old_hist.yrs <- hist.yrs[which(!(hist.yrs %in% new_hist.yrs))]
  
  
#  if(is.null(fbar)) fbar <- c(ages[1], ages[length(ages)])
  
  yrs <- hist.yrs[1]:sim.yrs[length(sim.yrs)]
  
  nages_stk  <- sapply(ages, function(x) length(x))
  nyear <- length(yrs)
  
  # Fcube format catch and effort (it can be excel, cvs or RData, but both in the same format)
  fmt <- strsplit(catch_obj,'.', fixed = TRUE)[[1]][length(strsplit(catch_obj,'.', fixed = TRUE)[[1]])] 
  if(fmt %in% c('xls', 'xlsx')){
    catch_wb  <- loadWorkbook(catch_obj)
    effort_wb <- loadWorkbook(effort_obj)
    catch  <- readWorksheet(catch_wb, sheet = 1)
    # check that all required columns are available
    if (!all(c("year","fleet","metier","stock","category","catch") %in% names(catch)))
      stop(paste("Columns 'year', 'fleet', 'metier', 'stock', 'category' and 'catch' are required in ", catch_obj, " file",sep=''))
    # transform to character (if necessary)
    catch$fleet    <- as.character(catch$fleet)
    catch$metier   <- as.character(catch$metier)
    catch$stock    <- as.character(catch$stock)
    catch$category <- as.character(catch$category)
    effort <- readWorksheet(effort_wb, sheet = 1)
    # check that all required columns are available
    if (!all(c("year","fleet","metier","effort") %in% names(effort)))
      stop(paste("Columns 'year', 'fleet', 'metier' and 'effort' are required in ", effort_obj, " file",sep=''))
    # transform to character (if necessary)
    effort$fleet  <- as.character(effort$fleet)
    effort$metier <- as.character(effort$metier)
  }
  else{
    if(fmt == 'RData'){
      load(catch_obj)
      load(effort_obj)
    }
    else{
      if(fmt == 'csv'){
        catch <- read.csv(catch_obj)
        effort <- read.csv(effort_obj)
      }
      else{
        stop('The format of catch and effort data must one of "RData", "csv", "xls" or "xlsx"')
      }
    }
  }

  # Subset only to the years in new_hist.yrs
  catch  <- subset(catch, year %in% new_hist.yrs)
  effort <- subset(effort, year %in% new_hist.yrs)
  
  # calculate total catch by stock and year
  catch_yr <- aggregate(catch ~ year + stock, catch, sum, na.rm = TRUE)
  
  # calculate total effort by fleet and year
  effort_yr <- aggregate(effort ~ year + fleet, effort, sum, na.rm = TRUE)
  
  # Aggregate by stock, fleet, metier and category just in case there are duplicates
  catch  <- aggregate(catch ~ stock + category + year + fleet + metier, catch, sum)
  effort <- aggregate(effort ~ year + fleet + metier, effort, sum)
  # add new column to the data.frame 'prop_mt' and 'prop_flmt' with the proportion by metier and by fleet and metier, respectively
  catch <- ddply(catch, c("stock", "category", "year", "fleet"), transform, prop_mt = catch/sum(catch))
  catch <- ddply(catch, c("stock", "category", "year"), transform, prop_flmt = catch/sum(catch))
  # when catches==0, prop=NaN  --> set to 0
  catch$prop_mt[is.nan(catch$prop_mt)]     <- 0
  catch$prop_flmt[is.nan(catch$prop_flmt)] <- 0
  # add a column to the data.frame 'prop' with the proportion by metier
  effort <- ddply(effort, c("year", "fleet"), transform, prop = effort/sum(effort))
  
  if(excel == TRUE){
    nit <- 1
  } else {
    # stk_data <- stk_obj
    flt_data <- flt_obj
    nit <- dim(data$n)[3]
    
    flt_nms        <- names(flt_data)
    flt_mt_nms     <- lapply(flt_data, function(x) names(x@metiers)) 
    flt_mt_stk_nms <- lapply(flt_data, function(x) lapply(x@metiers, function(y) names(y@catches)))

        stk_nms <- unique(unlist(flt_mt_stk_nms))
    ages <- list()
    i <- 1
    for(fl in names(flt_data))
      for(mt in names(flt_data[[fl]]@metiers))
        if(stk_nms[i] %in% names(flt_data[[fl]]@metiers[[mt]]@catches)) {
          print(paste(stk_nms[i], fl, mt, sep=' - '))
          ages[[stk_nms[i]]] <- dimnames(flt_data[[fl]]@metiers[[mt]]@catches[[stk_nms[i]]]@landings.n)[[1]]
          i <- i+1
          if (i>length(stk_nms)) {
            break()
          } else next()
        }
  }
  
  
  #---------------------------------------------------------------------
  ## Create the structure of FLFleets if it does not exist
  #---------------------------------------------------------------------
  
  if(is.null(flt_obj)){
    
    flfleets <- list() 
      
    eff_flq <- FLQuant(dim = c(1, length(yrs), 1,1,1,nit), dimnames = list(quant = 'all', year = yrs, iter = 1:nit))
    eff1_flq <- FLQuant(dim = c(1, length(yrs), 1,1,1,nit), dimnames = list(quant = 'all', year = yrs, iter = 1:nit))
    eff0_flq <- FLQuant(dim = c(1, length(yrs), 1,1,1,nit), dimnames = list(quant = 'all', year = yrs, iter = 1:nit))
    
    stk_flq <- FLQuant(dim = c(1, length(yrs), 1,1,1,nit), dimnames = list(age = 'all', year = yrs, iter = 1:nit))
    stk1_flq <- FLQuant(dim = c(1, length(yrs), 1,1,1,nit), dimnames = list(age = 'all', year = yrs, iter = 1:nit))
    stk0_flq <- FLQuant(dim = c(1, length(yrs), 1,1,1,nit), dimnames = list(age = 'all', year = yrs, iter = 1:nit))
    
    flcs <- vector('list', length(stk_objs))
    names(flcs) <- names(stk_nms)
    for(st in stk_nms){
      stka_flq <- FLQuant(dim = c(nages_stk[[st]], length(yrs), 1,1,1,nit), dimnames = list(age = ages_stk[[st]], year = yrs, iter = 1:nit))
      stka1_flq <- FLQuant(1, dim = c(nages_stk[[st]], length(yrs), 1,1,1,nit), dimnames = list(age = ages_stk[[st]], year = yrs, iter = 1:nit))
      stka0_flq <- FLQuant(0, dim = c(nages_stk[[st]], length(yrs), 1,1,1,nit), dimnames = list(age = ages_stk[[st]], year = yrs, iter = 1:nit))
  
      flcs[[st]] <- FLCatchExt(name = st, 
                desc = paste('data imported from', stk_nms[st]), 
                range = c(min = as.numeric(ages[[st]][1]), max = as.numeric(ages[[st]][length(ages[[st]])]), plusgroup =  as.numeric(ages[[st]][length(ages[[st]])]),  
                          minyear = as.numeric(hist.yrs[1]), maxyear = as.numeric(hist.yrs[length(hist.yrs)])),
                beta = stka1_flq,  
                alpha = stka1_flq,
                discards.n = stka0_flq, discards.sel = stka0_flq, discards = stk0_flq,
                landings.sel = stka1_flq)
    }
    
    for(fl in flt_nms){
      
      print(fl)
      
      flmt <-list()
      for(mt in flt_mt_nms[[fl]]){
        flmt[[mt]] <- FLMetierExt(name = mt, desc = 'data imported from', effshare = eff1_flq/length(flt_mt_nms[[fl]]), catches = FLCatchesExt(flcs[flt_mt_stk_nms[[fl]][[mt]]]))
        
      }
      flfleets[[fl]] <- FLFleetExt(name = fl, desc = 'data imported from', effort = eff1_flq, fcost = eff_flq, capacity = eff_flq, crewshare = eff0_flq,
              metiers = FLMetiersExt(flmt))
    }
  } else flfleets <- flt_obj
  
  
  #---------------------------------------------------------------------
  ## Checking files and content
  #---------------------------------------------------------------------
  
  # check that all files are available in the assigned paths
  if (sum(!caa_objs %in% dir(caa_objs_path))>0) 
    stop(paste("Following files missing in '", caa_objs_path, "': \n", paste(caa_objs[!(caa_objs %in% dir(caa_objs_path))], collapse = '\n'), sep=''))
  if(update_price == TRUE) {if (sum(!price_objs %in% dir(price_objs_path))>0) 
    stop(paste("Following files missing in '", price_objs_path,"': \n", paste(price_objs[!(price_objs %in% dir(price_objs_path))], collapse = '\n'), sep=''))}
  if (!file.exists( catch_obj)) stop("'catch_obj' file is not available")
  if (!file.exists( effort_obj)) stop("'effort_obj' file is not available")
  
  for (st in stk_nms) {
    # caa
    if(caaOpt[st] %in% c(1,2)){ 
      flst <- lapply(flt_mt_stk_nms, function(x) unique(unlist(x)))
      flst_caa <- flst
      for (i in 1:length(flst))
        flst_caa[[i]]   <- paste('caa_', names(flst)[i], "_", flst[[i]], ".xlsx", sep = "")
      files_caa   <- sort(unlist(flst_caa))
      if (sum(!files_caa %in% caa_objs)>0) 
        stop(paste("File names missing in 'caa_objs': \n", paste(files_caa[!files_caa %in% caa_objs], collapse = '\n'), sep=''))
    } else if (caaOpt[st] %in% c(3,4,5)) {
      file_caa <- paste( caa_objs_path, paste('caa_', st, ".xlsx", sep = ""), sep='/') 
      if (!file.exists(file_caa)) stop(paste("'", file_caa, "' file is not available", sep=''))
    }
    # prices
    if (update_price==TRUE) {
      if(priceOpt[st] %in% c(1,2)){ 
        flst <- lapply(flt_mt_stk_nms, function(x) unique(unlist(x)))
        flst_price <- flst
        for (i in 1:length(flst))
          flst_price[[i]] <- paste('price_', names(flst)[i], "_", flst[[i]], ".xlsx", sep = "")
        files_price <- sort(unlist(flst_price))
        if (sum(!files_price %in% price_objs)>0) 
          stop(paste("File names missing in 'price_objs': \n", paste(files_price[!files_price %in% price_objs], collapse = '\n'), sep=''))
      } else if (priceOpt[st] %in% c(3,4,5)) {
        file_price <- paste( price_objs_path, paste('price_', st, ".xlsx", sep = ""), sep='/')
        if (!file.exists(file_price)) stop(paste("'", file_price, "' file is not available", sep=''))
      }
    }
    
  }
  
  # Check that catch and effort files (catch_obj and effort_obj, respectively) 
  # contain the same fleets and metiers
  if ( sum(sort(flt_nms) != sort(unique(catch$fleet)))>0 )
    stop(paste("Check '", catch_obj, "' file, as only following fleets should appear: ", paste(flt_nms, collapse = ', '), sep = ''))
  if ( sum(sort(flt_nms) != sort(unique(effort$fleet)))>0 )
    stop(paste("Check '", effort_obj, "' file, as only following fleets should appear: ", paste(flt_nms, collapse = ', '), sep = ''))
  
  flmt_nms <- unique(catch[,c('fleet','metier')])
  
  for (fl in flt_nms) {
    
    if (sum(sort(flt_mt_nms[[fl]]) != sort(unique(subset(catch, fleet==fl)$metier)))>0)
      stop(paste("Check '", catch_obj, "' file, as for fleet '", fl, "' only following metiers should appear: ", paste(flt_mt_nms[[fl]], collapse = ', '), sep = ''))
    if (sum(sort(flt_mt_nms[[fl]]) != sort(unique(subset(effort, fleet==fl)$metier)))>0)
      stop(paste("Check '", effort_obj, "' file, as for fleet '", fl, "' only following metiers should appear: ", paste(flt_mt_nms[[fl]], collapse = ', '), sep = ''))
    
    for (mt in flt_mt_nms[[fl]]) {
    
      if (sum(sort(flt_mt_stk_nms[[fl]][[mt]]) != sort(unique(subset(catch, fleet==fl & metier==mt)$stock)))>0)
        stop(paste("Check '", catch_obj, "' file, as for fleet '", fl, "' and metier '", mt, "' only following stocks should appear: ", paste(flt_mt_stk_nms[[fl]][[mt]], collapse = ', '), sep = ''))
      if (sum(sort(flt_mt_stk_nms[[fl]][[mt]]) != sort(unique(subset(effort, fleet==fl & metier==mt)$stock)))>0)
        stop(paste("Check '", effort_obj, "' file, as for fleet '", fl, "' and metier '", mt, "' only following stocks should appear: ", paste(flt_mt_stk_nms[[fl]][[mt]], collapse = ', '), sep = ''))
      
    } # end mt
  } # end fl


  #---------------------------------------------------------------------
  ## Now, we have the structure and we need to fill in with the data.
  #---------------------------------------------------------------------
  ## Read in the weights at age
  #-----------------------------
  cat('--------------------------------------------------------------------\n')   
  cat(' Weights at age \n') 
  cat('--------------------------------------------------------------------\n') 
  
  wts.land <- wts.disc <- list()
  
  for(st in stk_nms){
    
   # browser()
    # print(st)
    fmt_stk <- strsplit(stk_objs[[st]], '.', fixed = TRUE)[[1]][length(strsplit(stk_objs[[st]], '.', fixed = TRUE)[[1]])]
    print(st)
    if(fmt_stk %in% c('xls', 'xlsx')){
      wb <- loadWorkbook(stk_objs[[st]], create = FALSE)
      sheets <- getSheets(wb)
      wl <- ifelse('wl' %in% sheets, 'wl', 'wt')
      wd <- ifelse('wd' %in% sheets, 'wd', wl)
      
      wts.land[[st]] <- as.matrix(readWorksheet(wb, sheet = wl, header = TRUE, startRow = 1, startCol = 2, 
                                                endRow = nages_stk[[st]] + 1))
      wts.disc[[st]] <- as.matrix(readWorksheet(wb, sheet = wd, header = TRUE, startRow = 1, startCol = 2, 
                                                endRow = nages_stk[[st]]  + 1))
      yrs_nms <- paste('X', hist.yrs, sep = "")
    }
    else{
      if(fmt_stk == 'RData'){
        
        loadToEnv(stk_objs[[st]])
        
        wl <- ifelse('wl' %in% names(data), 'wl', 'wt')
        wd <- ifelse('wd' %in% names(data), 'wd', wl)
        
        wts.land[[st]] <- apply(data[[wl]],1:2, median)
        wts.disc[[st]] <- apply(data[[wd]],1:2, median)
#if(st == 'MEG') browser()        
        yrs_nms <- colnames(data[[wl]])[colnames(data[[wl]]) %in% hist.yrs]
      }
      else{ stop('stock data must be provided in excel or R format')}
      }
    
    for(fl in names(flfleets)){
      for(mt in names(flfleets[[fl]]@metiers)){
        if(!(st %in% catchNames(flfleets[[fl]][[mt]]))) next

        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt[, hist.yrs] <- wts.land[[st]][, yrs_nms] 
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt[, hist.yrs] <- wts.disc[[st]][, yrs_nms] 
      }
    }
    
  }
  
  
  
  ## Landings and Discards at age
  #----------------------------------------------------------------
  cat('--------------------------------------------------------------------\n')   
  cat('All to do with catches \n') 
  cat('--------------------------------------------------------------------\n') 
  
 # browser()
  for(fl in names(flfleets)){
    
    for(mt in names(flfleets[[fl]]@metiers)){
      for(st in names(flfleets[[fl]]@metiers[[mt]]@catches)){
        xlcFreeMemory()
        fltmt <- paste(fl, mt, sep = "_")
        # print(fltmt)
        # print(st)
        
        # The CAA is provided by metier, for all the metiers
        #---------------------------------------------------
          if(caaOpt[st] == 1){ 

            caa_obj    <- paste("caa_", fl, "_", st, ".xlsx", sep = "")
            wb_caa     <- loadWorkbook(file.path(caa_objs_path, caa_obj))
            sheets_caa <- getSheets(wb_caa)
            
            if(!mt %in% sheets_caa)
              stop(paste("Sheet '", mt, "' missing in file: '", caa_obj, "'", sep=''))
            
            la <- as.matrix(readWorksheet(wb_caa, sheet = mt, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[st] + 1))
            da <- as.matrix(readWorksheet(wb_caa, sheet = mt, header = TRUE, startRow = nages_stk[st] +  3, startCol = 2, endRow = 2*nages_stk[st] + 3))
            
            la.yrs <-  ac(readWorksheet(wb_caa, sheet = mt, header = FALSE, startRow = 1, startCol = 2, endRow = 1))
            da.yrs <-  ac(readWorksheet(wb_caa, sheet = mt, header = FALSE, startRow = nages_stk[st] +  3, startCol = 2, endRow = nages_stk[st] +  3))
            
            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n[,la.yrs] <- la
            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n[,da.yrs] <- da
            
          }
          
        # The CAA is provided at fleet level.
        #------------------------------------
          if(caaOpt[st] == 2){
            
            caa_obj    <- paste("caa_", fl, "_", st, ".xlsx", sep = "")
            wb_caa     <- loadWorkbook(file.path(caa_objs_path, caa_obj))
            sheets_caa <- getSheets(wb_caa)
            
            if(!fl %in% sheets_caa)
              stop(paste("Sheet '", fl, "' missing in file: '", caa_obj, "'", sep=''))
            
            land_prop_mt <- subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$prop_mt
            disc_prop_mt <- subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$prop_mt
            
            land_prop_mt_yrs <- as.character(subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$year)
            disc_prop_mt_yrs <- as.character(subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$year)
            
            la <- as.matrix(readWorksheet(wb_caa, sheet = fl, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[[st]] + 1))
            da <- as.matrix(readWorksheet(wb_caa, sheet = fl, header = TRUE, startRow = nages_stk[[st]] + 3, startCol = 2, endRow = 2*nages_stk[[st]] + 3))
            
            la.yrs <-  readWorksheet(wb_caa, sheet = fl, header = FALSE, startRow = 1, startCol = 2, endRow = 1)
            da.yrs <-  readWorksheet(wb_caa, sheet = fl, header = FALSE, startRow = nages_stk[[st]] +  3, startCol = 2, endRow = nages_stk[[st]] +  3)
            
            colnames(la) <- la.yrs
            colnames(da) <- da.yrs

            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n[,land_prop_mt_yrs] <- sweep(la[,land_prop_mt_yrs, drop=FALSE],2,land_prop_mt,"*")
            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n[,disc_prop_mt_yrs] <- sweep(da[,disc_prop_mt_yrs, drop=FALSE],2,disc_prop_mt,"*")

          }
        
        # There is a correspondence between fleet segments in the CAA and fleet-metier disagregation.
        #--------------------------------------------------------------------------------------------
          if(caaOpt[st] == 3){
            
            wb_corres     <- loadWorkbook(caa_flt_mt_correspondences)
            sheets_corres <- getSheets(wb_corres)
            
            # identify the fleet_segment that corresponds for stock 'st' with fleet 'fl' and metier 'mt'.
            corres_st <-  readWorksheet(wb_corres, sheet = st, header = TRUE)
            
            caa_obj    <- paste("caa_", st, ".xlsx", sep = "")
            wb_caa     <- loadWorkbook(file.path(caa_objs_path,caa_obj))
            sheets_caa <- getSheets(wb_caa)
            
            fleetSeg  <- subset(corres_st, fleet_flbeia == fl & metier_flbeia == mt)[,3]
            
            if (identical(fleetSeg, character(0))) 
              stop(paste("Check '", caa_flt_mt_correspondences, "' file, as 'fleet_segment' correspondence missing for fleet '", fl, "' and metier '", mt, "'", sep=''))
            
            if(!fleetSeg %in% sheets_caa)
              stop(paste("Sheet '", fleetSeg, "' missing in file: '", caa_obj, "'", sep=''))
            
            land_mt <- subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$catch
            disc_mt <- subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$catch
            
            land_mt_yrs <- as.character(subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$year)
            disc_mt_yrs <- as.character(subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$year)
            
            cat(fl, mt, st, '\n')
          #  browser()
            
            la <- as.matrix(readWorksheet(wb_caa, sheet = fleetSeg, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[st] + 1))
            da <- as.matrix(readWorksheet(wb_caa, sheet = fleetSeg, header = TRUE, startRow = nages_stk[st] +  3, startCol = 2, endRow = 2*nages_stk[st] + 3))
            
            dy <- colnames(la)
            hy <- paste('X', hist.yrs, sep="")
            selyrs <- dy[which(dy %in% hy)]
            
            la <- la[,selyrs]
            if(dim(da)[1]!= 0) da <- da[,selyrs]
              
            la.yrs <-  ac(readWorksheet(wb_caa, sheet = fleetSeg, header = FALSE, startRow = 1, startCol = 2, endRow = 1))
            da.yrs <-  ac(readWorksheet(wb_caa, sheet = fleetSeg, header = FALSE, startRow = nages_stk[st] +  3, startCol = 2, endRow = nages_stk[st] +  3))
       
            law <- flfleets[[fl]][[mt]][[st]]@landings.wt[, substr(colnames(la),2,5), drop=T]
            daw <- flfleets[[fl]][[mt]][[st]]@discards.wt[, substr(colnames(da),2,5), drop=T] 
            
        #  browser()
            pla <- sweep(la*law, 2, apply(la*law,2,sum), "/") # catch proportions by age
            if(dim(da)[1] != 0){ pda <- sweep(da*daw, 2, apply(da*daw,2,sum), "/")} # catch proportions by age
            else  {pda <- sweep(da, 2, apply(da,2,sum), "/")}
            
            
            colnames(pla) <- substr(selyrs, 2,5)
            if(dim(pda)[1] != 0) colnames(pda) <- substr(selyrs, 2,5)
            
            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n[,land_mt_yrs] <- sweep(pla[,land_mt_yrs, drop=FALSE],2,land_mt,"*")/flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt[,land_mt_yrs,drop=T]
            if(dim(pda)[1] != 0){ 
            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n[,disc_mt_yrs] <- sweep(pda[,disc_mt_yrs, drop=FALSE],2,disc_mt,"*")/flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt[,disc_mt_yrs,drop=T]}
            else{flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n[,disc_mt_yrs] <- 0}
            
          }
        
        # The CAA is given at stock level
        #--------------------------------------------------------------------------------------------
          if(caaOpt[st] == 4){
          #  if(st == 'HKE')  browser()
            
            caa_obj    <- paste("caa_", st, ".xlsx", sep = "")
            wb_caa     <- loadWorkbook(file.path(caa_objs_path, caa_obj))
            sheets_caa <- getSheets(wb_caa)
            
            if(!st %in% sheets_caa)
              stop(paste("Sheet '", st, "' missing in file: '", caa_obj, "'", sep=''))
            
            land_prop_flmt <- subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$prop_flmt
            disc_prop_flmt <- subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$prop_flmt

            land_prop_flmt_yrs <- as.character(subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$year)
            disc_prop_flmt_yrs <- as.character(subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$year)

            la <- as.matrix(readWorksheet(wb_caa, sheet = st, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[[st]] + 1))
            da <- as.matrix(readWorksheet(wb_caa, sheet = st, header = TRUE, startRow = nages_stk[[st]] + 3, startCol = 2, endRow = 2*nages_stk[[st]] + 3))

            la.yrs <-  readWorksheet(wb_caa, sheet = st, header = FALSE, startRow = 1, startCol = 2, endRow = 1)
            da.yrs <-  readWorksheet(wb_caa, sheet = st, header = FALSE, startRow = nages_stk[[st]] +  3, startCol = 2, endRow = nages_stk[[st]] +  3)

            colnames(la) <- la.yrs
            colnames(da) <- da.yrs

            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n[,land_prop_flmt_yrs] <- sweep(la[,land_prop_flmt_yrs, drop=FALSE],2,land_prop_flmt,"*")
            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n[,disc_prop_flmt_yrs] <- sweep(da[,disc_prop_flmt_yrs, drop=FALSE],2,disc_prop_flmt,"*")
            
          }
        # The CAA is given at stock level but there is CAA[fl,mt] available in some historical years in the FLFleets obj.
        #---------------------------------------------------------------------------------------------------------------
        if(caaOpt[st] == 5){
          
         # browser()
          
          if(is.null(flt_obj)) stop('Option 5 cannot be used if flt_obj is not provided!')
          
          land_flmt <- subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$catch
          disc_flmt <- subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$catch
          
          land_flmt_yrs <- as.character(subset(catch, stock == st & category == 'landings' & metier == mt & fleet == fl)$year)
          disc_flmt_yrs <- as.character(subset(catch, stock == st & category == 'discards' & metier == mt & fleet == fl)$year)
          
          law <- (flt_obj[[fl]][[mt]][[st]]@landings.n[,old_hist.yrs]*flt_obj[[fl]][[mt]][[st]]@landings.wt[,old_hist.yrs])
          daw <- (flt_obj[[fl]][[mt]][[st]]@landings.n[,old_hist.yrs]*flt_obj[[fl]][[mt]][[st]]@landings.wt[,old_hist.yrs])
          
          plaa <- yearMeans(law%/%quantSums(law))
          pdaa <- yearMeans(daw%/%quantSums(daw))
          
          flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n[,new_hist.yrs] <- quantSums(plaa%*%FLQuant(land_flmt, dim = c(1,length(new_hist.yrs))))
          flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n[,new_hist.yrs] <- quantSums(pdaa%*%FLQuant(disc_flmt, dim = c(1,length(new_hist.yrs))))
          
        }

      }
    }
  }


  if(caaOpt[st] == 5){ # Correct the landings.n and discards.n of the fleets to equal the CAA in the stk_obj and in the flt_obj.
  
    if(is.null(flt_obj)) stop('Option 5 cannot be used if flt_obj is not provided!')
  
    wb_caa <- loadWorkbook(file.path(caa_objs_path, paste('caa_', st, ".xlsx", sep = "")))
    la <- as.matrix(readWorksheet(wb_caa, sheet = st, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[[st]] + 1))
    da <- as.matrix(readWorksheet(wb_caa, sheet = st, header = TRUE, startRow = nages_stk[[st]] + 3, startCol = 2, endRow = 2*nages_stk[[st]] + 3))
    la.yrs <-  readWorksheet(wb_caa, sheet = st, header = FALSE, startRow = 1, startCol = 2, endRow = 1)
    da.yrs <-  readWorksheet(wb_caa, sheet = st, header = FALSE, startRow = nages_stk[[st]] +  3, startCol = 2, endRow = nages_stk[[st]] +  3)
    colnames(la) <- la.yrs
    colnames(da) <- da.yrs
}
  
 
#  browser()
  
    ## Read in the prices at age
  #-----------------------------  
  cat('--------------------------------------------------------------------\n')   
  cat(' Prices \n') 
  cat('--------------------------------------------------------------------\n') 
if(update_price == TRUE){
  for(fl in names(flfleets)){
    
    for(mt in names(flfleets[[fl]]@metiers)){
      for(st in names(flfleets[[fl]]@metiers[[mt]]@catches)){
        xlcFreeMemory()
        fltmt <- paste(fl, mt, sep = "_")
        print(fltmt)
        print(st)

        # The price is provided by metier, for all the metiers
        if(priceOpt[st] == 1){ 
          
          price_obj    <- paste("price_", fl, "_", st, ".xlsx", sep = "")
          wb_price     <- loadWorkbook(file.path(price_objs_path, price_obj))
          sheets_price <- getSheets(wb_price)

          if(!mt %in% sheets_price)
            stop(paste("Sheet '", mt, "' missing in file: '", price_obj, "'", sep=''))
          
          pa <- as.matrix(readWorksheet(wb_price, sheet = mt, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[st] + 1))
         
          pa.yrs <-  ac(readWorksheet(wb_price, sheet = mt, header = FALSE, startRow = 1, startCol = 2, endRow = 1))

          flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price[,pa.yrs] <- pa

        }
        
        # The price is provided at fleet level.
        #------------------------------------
        if(priceOpt[st] == 2){
          
          # browser()
          
          price_obj    <- paste("price_", fl, "_", st, ".xlsx", sep = "")
          wb_price     <- loadWorkbook(file.path(price_objs_path, price_obj))
          sheets_price <- getSheets(wb_price)
          
          if(!fl %in% sheets_price)
            stop(paste("Sheet '", fl, "' missing in file: '", price_obj, "'", sep=''))
          
          pa <- as.matrix(readWorksheet(wb_price, sheet = fl, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[[st]] + 1))
          
          pa.yrs <-  ac(readWorksheet(wb_price, sheet = fl, header = FALSE, startRow = 1, startCol = 2, endRow = 1))
          
          flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price[,pa.yrs] <- pa

        }
        
        # There is a correspondence between fleet segments in the price and fleet-metier disagregation.
        #--------------------------------------------------------------------------------------------
        if(priceOpt[st] == 3){
          
          wb_corres     <- loadWorkbook(paa_flt_mt_correspondences)
          sheets_corres <- getSheets(wb_corres)
          
          # identify the fleet_segment that corresponds for stock 'st' with fleet 'fl' and metier 'mt'.
          corres_st <-  readWorksheet(wb_corres, sheet = st, header = TRUE)
          
          price_obj    <- paste("price_", st, ".xlsx", sep = "")
          wb_price     <- loadWorkbook(file.path(price_objs_path, price_obj))
          sheets_price <- getSheets(wb_price)
          
          fleetSeg  <- subset(corres_st, fleet_flbeia == fl & metier_flbeia == mt)[,3]
          
          if (identical(fleetSeg, character(0))) 
            stop(paste("Check '", caa_flt_mt_correspondences, "' file, as 'fleet_segment' correspondence missing for fleet '", fl, "' and metier '", mt, "'", sep=''))
          
          if(!fleetSeg %in% sheets_price)
            stop(paste("Sheet '", fleetSeg, "' missing in file: '", price_obj, "'", sep=''))
          
          pa <- as.matrix(readWorksheet(wb_price, sheet = fleetSeg, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[st] + 1))
          
          pa.yrs <-  ac(readWorksheet(wb_price, sheet = fleetSeg, header = FALSE, startRow = 1, startCol = 2, endRow = 1))
          
          flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price[,pa.yrs] <- pa
          
        }
        
        # The price is given at fleet level
        #--------------------------------------------------------------------------------------------
        if(priceOpt[st] == 4){

          price_obj    <- paste("price_", st, ".xlsx", sep = "")
          wb_price <- loadWorkbook(file.path(price_objs_path, price_obj))
          sheets_price <- getSheets(wb_price)
          
          if(!st %in% sheets_price)
            stop(paste("Sheet '", st, "' missing in file: '", price_obj, "'", sep=''))
          
          pa <- as.matrix(readWorksheet(wb_price, sheet = st, header = TRUE, startRow = 1, startCol = 2, endRow = nages_stk[[st]] + 1))
          
          pa.yrs <-  ac(readWorksheet(wb_price, sheet = st, header = FALSE, startRow = 1, startCol = 2, endRow = 1))
          
          flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price[,pa.yrs] <- pa
          
        }

      }
    }
  }
  
}  
  # browser()
  
  ## Calculate and fill in Effort share, vcost & effort
  #-------------------------------------------
  cat('--------------------------------------------------------------------\n')   
  cat(' Effort, effort share and vcost \n') 
  cat('--------------------------------------------------------------------\n') 
#browser()
  for(fl in names(flfleets)){
    
    eff_fl     <- subset(effort, fleet == fl)
    eff.yrs <- as.character(unique(eff_fl$year))
    
    flfleets[[fl]]@effort[,eff.yrs] <- aggregate(effort ~ year, eff_fl, sum)$effort
    
    for(mt in names(flfleets[[fl]]@metiers)){
      cat(fl, mt, '\n')
      eff_mt <- subset(effort, fleet == fl & metier == mt)
      flfleets[[fl]]@metiers[[mt]]@effshare[,as.character(eff_mt$year)] <- eff_mt$prop
    }
    
  }
  

  ## Fill  fcost, capacity & crewshare
  #-----------------------------------------------
  
  
  ## Values for proyection
  #-----------------------------------------------
  
  cat('--------------------------------------------------------------------\n')   
  cat(' Set values for projection period \n') 
  cat('--------------------------------------------------------------------\n') 

  for(fl in names(flfleets)){
    
    cat('---------------------', fl,'fleet,','------------------\n')
    
    # effort
    effort(flfleets[[fl]])[,sim.yrs,]   <-  yearMeans(effort(flfleets[[fl]])[,mean.yrs,])
    
    if(any(is.na(effort(flfleets[[fl]])[,mean.yrs,]))) {
      cat(paste("warning: NAs in effort for 'mean.yrs' and fleet '", fl, "' \n", sep = ''))
      if(any(is.na(effort(flfleets[[fl]])[,sim.yrs,]))) 
        cat(paste("warning: all NAs in effort for 'sim.yrs' and fleet '", fl, "' \n", sep = ''))
    }
    
    # fcost
    flfleets[[fl]]@fcost[,sim.yrs,]     <-  yearMeans(flfleets[[fl]]@fcost[,mean.yrs,])
    if(any(is.na(flfleets[[fl]]@fcost[,mean.yrs,]))) {
      cat(paste("warning: NAs in fcost for 'mean.yrs' and fleet '", fl, "' \n", sep = ''))
      if (any(is.na(flfleets[[fl]]@fcost[,sim.yrs,]))) 
        cat(paste("warning: all NAs in fcost for 'sim.yrs' and fleet '", fl, "' \n", sep = ''))
    }
    
    # capacity
    flfleets[[fl]]@capacity[,sim.yrs,]  <-  yearMeans(flfleets[[fl]]@capacity[,mean.yrs,])
    if(any(is.na(flfleets[[fl]]@capacity[,mean.yrs,]))) {
      cat(paste("warning: NAs in capacity for 'mean.yrs' and fleet '", fl, "' \n", sep = ''))
      if (any(is.na(flfleets[[fl]]@capacity[,sim.yrs,])))
        cat(paste("warning: all NAs in capacity for 'sim.yrs' and fleet '", fl, "' \n", sep = ''))
    }
    
    # crewshare
    flfleets[[fl]]@crewshare[,sim.yrs,] <-  yearMeans(flfleets[[fl]]@crewshare[,mean.yrs,])
    if(any(is.na(flfleets[[fl]]@crewshare[,mean.yrs,]))) {
      cat(paste("warning: NAs in crewshare for 'mean.yrs' and fleet '", fl, "' \n", sep = ''))
      if (any(is.na(flfleets[[fl]]@crewshare[,sim.yrs,]))) 
        cat(paste("warning: all NAs in crewshare for 'sim.yrs' and fleet '", fl, "' \n", sep = ''))
    }
    
    all.efs <- 0
    for(mt in names(flfleets[[fl]]@metiers)){
      
      cat('---------------------', fl,'fleet,',mt,' metier,','---------------------\n')
      
      # effshare
      met.efs <- yearMeans(flfleets[[fl]]@metiers[[mt]]@effshare[,mean.yrs,])
      if (fl==names(flfleets)[length(names(flfleets))] & length(names(flfleets))>1){ 
        # the sum of all effshare must be one
        flfleets[[fl]]@metiers[[mt]]@effshare[,sim.yrs,] <- 1 - all.efs
      } else {
        flfleets[[fl]]@metiers[[mt]]@effshare[,sim.yrs,] <-  met.efs
      } 
      all.efs <- all.efs + met.efs

      if(any(is.na(flfleets[[fl]]@metiers[[mt]]@effshare)[,mean.yrs,,,,])) {
        cat(paste("warning: NAs in effshare for 'mean.yrs', fleet '", fl, "' and metier '", mt, "' \n", sep = ''))
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@effshare)[,sim.yrs,,,,]))
          cat(paste("warning: all NAs in effshare for 'sim.yrs', fleet '", fl, "' and metier '", mt, "' \n", sep = ''))
      }

      # vcost
      flfleets[[fl]]@metiers[[mt]]@vcost[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@vcost[,mean.yrs,])
      if(any(is.na(flfleets[[fl]]@metiers[[mt]]@vcost)[,mean.yrs,,,,])) {
        cat(paste("warning: NAs in vcost for 'mean.yrs', fleet '", fl, "' and metier '", mt, "' \n", sep = ''))
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@vcost)[,sim.yrs,,,,]))
          cat(paste("warning: all NAs in vcost for 'sim.yrs', fleet '", fl, "' and metier '", mt, "' \n", sep = ''))
      }
#browser()
      for(st in names(flfleets[[fl]]@metiers[[mt]]@catches)){
        
        cat('---------------------', fl,'fleet,',mt,' metier,',st,'stock','---------------------\n')
        
        # landings.wt
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt[,mean.yrs,])
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt)[,mean.yrs,,,,])) {
          cat(paste("warning: NAs in landings.wt for 'mean.yrs', fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
          if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt)[,sim.yrs,,,,]))
            cat(paste("warning: NAs in landings.wt for 'sim.yrs', fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        }
        
        # landings
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings <- quantSums( flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.n*
                                                                            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.wt)
        
        # landings.sel
        if (any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel)[,mean.yrs,,,,])) {
          cat(paste("warning: NAs in landings.sel for fleet '", fl, "', metier '", mt, "' and stock '", st,".\n", sep = ''))
                  #  "', these have been replaced by 1 for computing means. \n", sep = ''))
         # flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,mean.yrs,][is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,mean.yrs,])] <- 1
        }
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel[,mean.yrs,])
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@landings.sel)[,sim.yrs,,,,]))
          cat(paste("warning: NAs in landings.sel projection for fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))

        # discards.wt
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt[,mean.yrs,])
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt)[,mean.yrs,,,,])) {
          cat(paste("warning: NAs in discards.wt for 'mean.yrs', fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
          if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt)[,sim.yrs,,,,]))
            cat(paste("warning: NAs in discards.wt for 'sim.yrs', fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        }
        
        # discards
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards <- quantSums( flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.n*
                                                                            flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.wt)
        
        # discards.sel
        if (any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel)[,mean.yrs,,,,])) {
          cat(paste("warning: NAs in discards.sel for fleet '", fl, "', metier '", mt, "' and stock '", st,"\n", sep = ''))
                   # "', these have been replaced by 0 for computing means. \n", sep = ''))
         # flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel[,mean.yrs,][is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel[,mean.yrs,])] <- 0
        }
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel[,mean.yrs,])
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@discards.sel)[,sim.yrs,,,,]))
          cat(paste("warning: NAs in discards.sel projection for fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        
        # alpha
        # if(any(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha[,mean.yrs,,,,]<0, na.rm = TRUE)) {
        #   stop(paste("Negative values in alpha for fleet '", fl, "', metier '", mt, "' and stock '", st, 
        #              "' in years used for computing means (i.e. 'mean.yrs'). \n", sep = ''))
        # } else if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha)[,mean.yrs,,,,])) {
        #   cat(paste("warning: NAs in alpha for fleet '", fl, "', metier '", mt, "' and stock '", st,
        #             "', these have been replaced by 1 for computing means. \n", sep = ''))
        #   flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha[,mean.yrs,][is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha[,mean.yrs,])] <- 1
        # }
        # flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha[,mean.yrs,])
        # if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@alpha)[,sim.yrs,,,,]))
        #   cat(paste("warning: NAs in alpha projection for fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        
        # beta
        # if(any(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta[,mean.yrs,,,,]<0, na.rm = TRUE)) {
        #   stop(paste("Negative values in beta for fleet '", fl, "', metier '", mt, "' and stock '", st, 
        #              "' in years used for computing means (i.e. 'mean.yrs'). \n", sep = ''))
        # } else if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta)[,mean.yrs,,,,])) {
        #   cat(paste("warning: NAs in beta for fleet '", fl, "', metier '", mt, "' and stock '", st,
        #             "', these have been replaced by 1 for computing means. \n", sep = ''))
        #   flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta[,mean.yrs,][is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta[,mean.yrs,])] <- 1
        # }
        # flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta[,mean.yrs,])
        # if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@beta)[,sim.yrs,,,,]))
        #   cat(paste("warning: NAs in beta projection for fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        # 
        
        # catch.q
        # if(any(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,mean.yrs,,,,]<0, na.rm = TRUE)) {
        #   stop(paste("Negative values in catch.q for fleet '", fl, "', metier '", mt, "' and stock '", st, 
        #              "' in years used for computing means (i.e. 'mean.yrs'). \n", sep = ''))
        # } else if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q)[,mean.yrs,,,,])) {
        #   cat(paste("warning: NAs in catch.q for fleet '", fl, "', metier '", mt, "' and stock '", st,
        #             "', these have been replaced by 0 for computing means. \n", sep = ''))
        #   flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,mean.yrs,][is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,mean.yrs,])] <- 0
        # }
        # flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q[,mean.yrs,])
        # if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q)[,sim.yrs,,,,]))
        #   cat(paste("warning: NAs in catch.q projection for fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        # 
        # price
        flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price[,sim.yrs,] <- yearMeans(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price[,mean.yrs,])
        if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price)[,mean.yrs,,,,])) {
          cat(paste("warning: NAs in price for 'mean.yrs', fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
          if(any(is.na(flfleets[[fl]]@metiers[[mt]]@catches[[st]]@price)[,sim.yrs,,,,]))
            cat(paste("warning: NAs in price for 'sim.yrs', fleet '", fl, "', metier '", mt, "' and stock '", st, "' \n", sep = ''))
        }
        
      } # END st
      
    } # END mt
    
    #! check sum(effshare by fleet) == 1
    # x <- lapply( flfleets[[fl]]@metiers, function(x) x@effshare)
    # Reduce('+',x)
    # Reduce('+',x)==0 | Reduce('+',x)==1
    
  } # END fl


  ## Return output
  #-----------------------------------------------
  
  flfleets <- FLFleetsExt(flfleets)
  
  return(flfleets)
         
}

