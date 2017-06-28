###############################################################################
# AUTHOR(DATE): Agurtzane Urtizberea, Dorleta Garcia and Sonia Sanchez
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        create.BDs.data
# NOTE #1:      Return a list of FLBDsim object called BDs
###############################################################################
#-------------------------------------------------------------------------
#' 
#' FLBEIA easy conditioning: BDs argument creator
#' 
#' create.BDs.data function creates a list of FLBDsim objects.
#' 
#' @param   ni Number of iterations (number).
#' @param   ns Number of seasons (number).
#' @param   yrs A vector with c(first.yr,proj.yr, last.yr) where:
#'  \itemize{
#'      \item first.yr: First year of simulation (number).
#'      \item proj.yr: First year of projection (number).
#'      \item last.yr: Last year of projection (number).}
#' @param   stks.data A list with the name of the stks and the following elements:
#'  \itemize{
#'      \item stk.unit: Number of units of the stock (number).
#'      \item stk.age.min: Minimum age class of the stock (number).
#'      \item stk.age.max: Maximum age class of the stock (number).
#'      \item stk_bd.model: Name of the model to simulate biomass dinamics of the stock (character).
#'      \item stk_params.name: Name of the parameters (vector).
#'      \item stk_params.array:	Parameter values (array).
#'      \item stk_biomass.flq: Biomass values (FLQuant).
#'      \item stk_catch.flq: Catch values (FLQuant).
#'      \item stk_range.plusgroup: Plusgroup age (numeric).
#'      \item stk_range.minyear: Minimum year (numeric).
#'      \item stk_alpha: Maximum variability of carrying capacity.}
#'  Optionals:
#'  \itemize{
#'      \item stk_gB.flq: Surplus production (FLQuant).
#'      \item stk_uncertainty.flq: Uncertainty (FLQuant).}
#'   
#' @return A list of FLBDsim objects.
#'     
#-------------------------------------------------------------------------------

create.BDs.data <- function (yrs,ns,ni,stks.data)
{
  ind <- unlist(sapply(stks.data,function(x) grep(x, pattern="_bd.model",value=TRUE)))
  nmstks <- unique(sub('.*?^(.*?)_bd.model*', '\\1', ind))
  BDs <- NULL
  
  if (length(nmstks) != 0) {
    n.stk.BD <- length(nmstks)
    
    first.yr <- yrs[["first.yr"]]
    proj.yr  <- yrs[["proj.yr"]]
    last.yr  <- yrs[["last.yr"]]
    proj.yrs       <- as.character(proj.yr:last.yr)
    hist.yrs       <- as.character(first.yr:(proj.yr-1))
    ny <- length(first.yr:last.yr)
    list.stks.unit <- lapply(stks.data, function(ch) grep(pattern="unit", ch, value = TRUE))
    list.stks.flq <- create.list.stks.flq(nmstks,yrs,ni,ns,list.stks.unit)  
    
    list.BDs <- list()
    for (i in 1:n.stk.BD) {
      nmstk <- nmstks[i]
      cat("=============", nmstk, "BD", "=============\n")
      flq.stk <- list.stks.flq[[nmstk]][, , 1]
      
      stk.model  <- get(grep(stks.data[[nmstk]],pattern="_bd.model", value = TRUE))
      stk.unit  <- get(grep(stks.data[[nmstk]],pattern=".unit", value = TRUE))
      stk.biomass  <- get(grep(stks.data[[nmstk]],pattern="_biomass.flq", value = TRUE))
      stk.catch  <- get(grep(stks.data[[nmstk]],pattern="_catch.flq", value = TRUE))
      stk.params.name  <- get(grep(stks.data[[nmstk]],pattern="_params.name", value = TRUE))
      stk.params  <- get(grep(stks.data[[nmstk]],pattern="_params.array", value = TRUE))
      stk.range.min       <- get(grep(stks.data[[nmstk]],pattern=".age.min", value = TRUE))
      stk.range.max       <- get(grep(stks.data[[nmstk]],pattern=".age.max", value = TRUE))
      stk.range.plusgroup       <- get(grep(stks.data[[nmstk]],pattern="_range.plusgroup", value = TRUE))
      stk.range.minyear       <- get(grep(stks.data[[nmstk]],pattern="_range.minyear", value = TRUE))
      stk.uncertainty       <- mget(grep(stks.data[[nmstk]],pattern="_uncertainty.flq", value = TRUE),envir=as.environment(1))
      stk.gB       <- mget(grep(stks.data[[nmstk]],pattern="_gB.flq", value = TRUE),envir=as.environment(1))
      
      if(length(stk.uncertainty)==0) stk.uncertainty  <- NA 
      if(length(stk.gB)==0) stk.gB  <- NA    
      
      stk.alpha      <- get(grep(stks.data[[nmstk]],pattern="_alpha", value = TRUE),envir=as.environment(1))


      params <- array(dim = c(length(stk.params.name), ny, ns, ni),
                dimnames = list(param = stk.params.name, year = ac(first.yr:last.yr),
                  season = ac(1:ns), iter = 1:ni))
      stk.bd <- FLBDsim(name = nmstk, model = stk.model,
                biomass = flq.stk, gB=flq.stk, catch = flq.stk, uncertainty = flq.stk,
                params = params)
      dimnames(stk.bd@params)$param <- stk.params.name
      
      stk.bd@alpha[] <- stk.alpha   
      stk.bd@range[["min"]] <- stk.range.min
      stk.bd@range[["max"]] <- stk.range.max
      stk.bd@range[["plusgroup"]] <- stk.range.plusgroup
      stk.bd@range[["minyear"]] <- stk.range.minyear
      stk.bd@range[["maxyear"]] <- last.yr
      if (!all(is.na(stk.uncertainty))) {
        stk.uncertainty <- stk.uncertainty[[1]]
        log.dim <- equal.flq.Dimnames(lflq = list(stk.uncertainty,
        stk.bd@uncertainty), 2)
      if (!log.dim)
           stop("BD uncertainty dimension names \n")
      if (!(any(dim(stk.uncertainty)[3] == c(1, stk.unit))))
            stop("in uncertainty number of stock units 1 or stk.unit")
      if (!(any(dim(stk.uncertainty)[4] == c(1, ns))))
            stop("in uncertainty number of seasons 1 or ns")
      if (!(any(dim(stk.uncertainty)[6] == c(1, ni))))
             stop("in uncertainty number of iterations 1 or ni")
            }
      else {
          stk.uncertainty = 1
          cat("BD uncertainty = 1 \n")
      }
      if (!all(is.na(stk.gB))) {
        stk.gB <- stk.gB[[1]]
        log.dim <- equal.flq.Dimnames(lflq = list(stk.gB,
                                                  stk.bd@gB[,hist.yrs]), 2)
        if (!log.dim)
          stop("BD gB dimension names \n")
        if (!(any(dim(stk.gB)[3] == c(1, stk.unit))))
          stop("in gB number of stock units 1 or stk.unit")
        if (!(any(dim(stk.gB)[4] == c(1, ns))))
          stop("in gB number of seasons 1 or ns")
        if (!(any(dim(stk.gB)[6] == c(1, ni))))
          stop("in gB number of iterations 1 or ni")
      }else {
        stk.gB = stk.bd@gB
        cat("gB is all NA-s")
      }
      if (!all(is.na(stk.biomass))) {
          log.dim <- equal.flq.Dimnames(lflq = list(stk.biomass,
              stk.bd@biomass[, hist.yrs]), 2)
          if (!log.dim)
              stop("in BD biomass dimension names \n")
          if (!(any(dim(stk.biomass)[3] == c(1, stk.unit))))
             stop("in biomass number of stock units 1 or stk.unit")
          if (!(any(dim(stk.biomass)[4] == c(1, ns))))
              stop("in biomass number of seasons 1 or ns")
          if (!(any(dim(stk.biomass)[6] == c(1, ni))))
              stop("in biomass number of iterations 1 or ni")
            }
       else {
           cat("BD biomass values all NA-s \n")
            }
           if (!all(is.na(stk.catch))) {
              log.dim <- equal.flq.Dimnames(lflq = list(stk.catch,
                  stk.bd@catch[, hist.yrs]), 2)
              if (!log.dim)
                  stop("in BD catch dimension names \n")
              if (!(any(dim(stk.catch)[3] == c(1, stk.unit))))
                  stop("in catch number of stock units 1 or stk.unit")
              if (!(any(dim(stk.catch)[4] == c(1, ns))))
                  stop("in catch number of seasons 1 or ns")
              if (!(any(dim(stk.catch)[6] == c(1, ni))))
                  stop("in catch number of iterations 1 or ni")
            }
            else {
                cat("BD catch values all NA-s \n")
            }
            log.dim <- equal.flq.Dimnames(lflq = list(stk.params,
                stk.bd@params), 1:4)
            if (!log.dim)
                stop("BD parameters dimension names \n")
            stk.bd@biomass[, hist.yrs] <- stk.biomass
            stk.bd@catch[, hist.yrs] <- stk.catch
            stk.bd@params <- stk.params
            dimnames(stk.bd@params)$param <- stk.params.name
            if(stk.bd@model=="PellaTom"){
              p <- stk.bd@params["p",,,]
              r <- stk.bd@params["r",,,]
              K <- stk.bd@params["K",,,]

              if(any(stk.bd@alpha<1) || any(stk.bd@alpha > min((p/r+1)^(1/p)))){
                stop("alpha<1 or alpha > min((p/r+1)^(1/p))")
              } 
            } 
            stk.bd@uncertainty[] <- stk.uncertainty
            stk.bd@gB[,hist.yrs] <- stk.gB
            
            if (!any(is.na(stk.params[, proj.yrs, , ]))) {
                if (!all(dim(stk.params) == dim(stk.bd@params))) {
                  stop("in BD parameters dimension names \n")
                }
            }
            else {
                stop("BD parameters all NA-s \n")
            }
            if (any(is.na(stk.bd@uncertainty[, proj.yrs]))) {
                stop("Na values in uncertainty in the projection years")
            }
            list.BDs[[i]] <- stk.bd
        }
        names(list.BDs) <- nmstks
        BDs <- list.BDs
    }
    return(BDs)
}
