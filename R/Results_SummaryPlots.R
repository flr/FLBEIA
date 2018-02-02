#------------------------------------------------------------------------------#
# Summary plots for output from summary functions
#
# (note input object should come from summary functions with quantiles)
# 
#------------------------------------------------------------------------------#
#' Summary plots of the FLBEIA output
#' 
#' Summarize the results in a plot.
#'
#' @details 
#' 
#'\itemize{
#'      \item{plotbioSum}: Plot summarising information on stock. 
#'         With one plot for each of the following indicators: "catch", "rec", "f" and "ssb".
#'         Input object should have the same format as the output of bioSum function.
#'      \item{plotfltSum}: Plot summarising information on fleet's economic indicators. 
#'         With one plot for each of the following indicators: "catch", "effort", "grossValue" and "grossSurplus".
#'         Input object should have the same format as the output of fltSum function.
#'}
#'
#' @param obj An object with the same format as bioSum and fltSum outputs, respectively.
#' @param stk.nam Character with the name of the stock for which summary information will be plotted. 
#'   If not defined, then the first one in obj will be selected by default.
#' @param flt.nam Character with the name of the fleet for which summary information will be plotted.
#'   If not defined, then the first one in obj will be selected by default.
#' @param Blim Numeric. Blim reference point for the stock (optional argument).=NA, Bpa=NA, proj.yr=NA) 
#' @param Bpa Numeric. Bpa reference point for the stock (optional argument).
#' @param proj.yr Numeric. Year in which projection period starts (optional argument). 
#' 
#' @return Plot. 
#'
#' @examples
#'\dontrun{
#'
#' library(FLBEIA)
#'
#' # Apply the summary plots to the examples runs in FLBEIA help page.
#' 
#' data(res_flbeia)
# 
#' #------------------------------------------------
#' # Example One: One stock, one fleet, one iter.
#' #------------------------------------------------
#' 
#' # Biological indicators.
#' plotbioSum( bioSum(oneRes), Blim=800, Bpa=1200, proj.yr=2010) # bioSum output in wide format
#' plotbioSum( bioSum(oneRes, long = FALSE))   # bioSum output in long format
#' plotbioSum( bioSum(oneRes), stk.nam='stk0') # if incorrect name for the stock
#' 
#' # Indicators at fleet level.
#' plotfltSum( fltSum(oneRes), proj.yr=2010) # fltSum output in wide format
#' plotfltSum( fltSum(oneRes, long = FALSE)) # fltSum output in long format
#' plotfltSum( fltSum(oneRes), flt.nam='stk1') # if incorrect name for the fleet
#' plotfltSum( fltSum(oneRes, byyear = FALSE)) # although seasonal disagregation, it is summarised by year
#' 
#' #------------------------------------------------
#' # Example OneIters: As one but with iterations.
#' #------------------------------------------------
#' 
#' # Biological indicators.
#' plotbioSum( bioSum(s1), Blim=800, Bpa=1200, proj.yr=2010) # bioSum output in wide format
#' plotbioSum( bioSum(s1, long = FALSE))   # bioSum output in long format
#' plotbioSum( bioSum(s1), stk.nam='stk0') # if incorrect name for the stock
#' 
#' 
#' # Indicators at fleet level.
#' plotfltSum( fltSum(s1), proj.yr=2010) # fltSum output in wide format
#' plotfltSum( fltSum(s1, long = FALSE)) # fltSum output in long format
#' plotfltSum( fltSum(s1), flt.nam='stk1') # if incorrect name for the fleet
#' plotfltSum( fltSum(s1, byyear = FALSE)) # although seasonal disagregation, it is summarised by year
#'
#' # also possible to plot information on various scenarios
#' sc11_bio <- bioSum(s1)
#' sc12_bio <- bioSum(s1, scenario='alt'); sc12_bio$value <- sc12_bio$value*1.05
#' plotbioSum(rbind(sc11_bio, sc12_bio), Blim=800, Bpa=1200, proj.yr=2010)
#' 
#' #------------------------------------------------
#' # Example Multi: Two stock, two fleet, four iters.
#' #------------------------------------------------
#' 
#' for (st in names(s2$stocks)) # one plot for each stock
#'   plotbioSum( bioSum(s2, scenario='s2'), stk.nam=st, proj.yr=2010)
#' 
#' for (fl in names(s2$fleets)) # one plot for each fleet
#'   plotfltSum( fltSum(s2, scenario='s2'), flt.nam=fl, proj.yr=2010)
#'   
#'}


#----------------------------------------------------------------------
# plotbioSum ( obj, stk.nam, Blim=NA, Bpa=NA, proj.yr=NA)
# 
# obj     = same format as output of bioSum function.
# stk.nam = name of the stock to be summarised.
# Bpa     = numeric (optional), with the precautionary biomass for the stock.
# Blim    = numeric (optional), with the limit biomass for the stock.
# proj.yr = numeric (optional), with the projection starting year.
#----------------------------------------------------------------------
#' @rdname plotbioSum
plotbioSum <- function( obj, stk.nam, Blim=NA, Bpa=NA, proj.yr=NA) {
  
  # Required indicators
  inds <- c('catch','rec','f','ssb')
  
  # Subset data for the selected stock
  if (missing(stk.nam)) { 
    stk.nam <- unique(obj$stock)[1]
  } else if (!stk.nam %in% unique(obj$stock)) {
    stop(paste("Stock '", stk.nam,"' does not exist. Options are: ", paste(unique(obj$stock), collapse = ", "), sep=''))
  }
  obj <- subset( obj, stock==stk.nam)
  
  # from wide to long format if necessary
  if (dim(obj)[2] > 6) {
    obj <- reshape( obj, direction = "long", varying=list(names(obj)[-c(1:4)]), 
                    v.names='value', idvar=names(obj)[c(1:4)], timevar='indicator', 
                    times=names(obj)[-c(1:4)])
  }
  
  # check the format of the object
  if(!"scenario" %in% colnames(obj)) {
    stop("Column 'scenario' is missing")
  } else if (!"year" %in% colnames(obj)) {
    stop("Column 'year' is missing")
  } else if (!"stock" %in% colnames(obj)) {
    stop("Column 'stock' is missing")
  } else if (!"iter" %in% colnames(obj)) {
    stop("Column 'iter' is missing")
  } else if (!"value" %in% colnames(obj)) {
    stop("Column 'value' is missing")
  } else if (sum(!inds %in% unique(obj$indicator))>0) {
    stop('Information on any of the required indicators (catch, rec, f, ssb) is missing')
  }
  
  # plotting
  d <- subset(bioSumQ(obj), indicator %in% inds)
  d$indicator <- factor( d$indicator, levels=inds)
  p <- ggplot( data=d, aes(x=year, y=q50, color=scenario)) + 
    geom_line() + 
    geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=scenario), alpha=0.5) + 
    facet_wrap(~ indicator, scales="free") + 
    theme_bw() + 
    theme(text=element_text(size=20),
          title=element_text(size=20,face="bold"),
          strip.text=element_text(size=20)) + 
    ylab("") + 
    ggtitle(stk.nam) + theme(plot.title = element_text(hjust = 0.5))
  if (!is.na(proj.yr))
    p <- p + geom_vline(xintercept = proj.yr, linetype = "longdash")
  if (!is.na(Blim) | !is.na(Bpa))
    p <- p + geom_hline(lty=2, data=data.frame(indicator=rep('ssb',2), Bref=c(Blim,Bpa)), 
                        aes(yintercept=Bref), color=c('red','orange'))
  print(p)
  
}


#----------------------------------------------------------------------
# plotfltSum ( obj, flt.nam, proj.yr=NA)
# 
# obj     = same format as output of bioSum function.
# flt.nam = name of the fleet to be summarised.
# proj.yr = numeric (optional), with the projection starting year.
#----------------------------------------------------------------------
#' @rdname plotbioSum
plotfltSum <- function( obj, flt.nam, proj.yr=NA) {
  
  # Required indicators
  inds <- c('catch','effort','grossValue','grossSurplus')
  
  # Subset data for the selected fleet
  if (missing(flt.nam)) {
    flt.nam <- unique(obj$fleet)[1]
  } else if (!flt.nam %in% unique(obj$fleet)){
    stop(paste("Fleet '", flt.nam,"' does not exist. Options are: ", paste(unique(obj$fleet), collapse = ", "), sep=''))
  }
  obj <- subset( obj, fleet==flt.nam)
  
  # from wide to long format if necessary
  if (dim(obj)[2] > 6 & !"season" %in% colnames(obj) | dim(obj)[2] > 7 & "season" %in% colnames(obj)) {
    obj <- reshape( obj, direction = "long", varying=list(names(obj)[-c(1:4)]), 
                    v.names='value', idvar=names(obj)[c(1:4)], timevar='indicator', 
                    times=names(obj)[-c(1:4)])
  }
  
  # check the format of the object
  if(!"scenario" %in% colnames(obj)) {
    stop("Column 'scenario' is missing")
  } else if (!"year" %in% colnames(obj)) {
    stop("Column 'year' is missing")
  } else if (!"fleet" %in% colnames(obj)) {
    stop("Column 'fleet' is missing")
  } else if (!"iter" %in% colnames(obj)) {
    stop("Column 'iter' is missing")
  } else if (!"value" %in% colnames(obj)) {
    stop("Column 'value' is missing")
  } else if (sum(!inds %in% unique(obj$indicator))>0) {
    stop('Information on any of the required indicators (catch, effort, grossValue, grossSurplus) is missing')
  }
  
  # plotting
  d <- subset(fltSumQ(obj), indicator %in% inds)
  d$indicator <- factor( d$indicator, levels=inds)
  d$year <- as.numeric(d$year)
  p <- ggplot( data=d, aes(x=year, y=q50, color=scenario)) + 
    geom_line() + 
    geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=scenario), alpha=0.5) + 
    facet_wrap(~ indicator, scales="free") + 
    # geom_vline(xintercept = oneMainC$sim.years[['initial']], linetype = "longdash") + 
    theme_bw() + 
    theme(text=element_text(size=20),
          title=element_text(size=20,face="bold"),
          strip.text=element_text(size=20)) + 
    ylab("") + 
    ggtitle(flt.nam) + theme(plot.title = element_text(hjust = 0.5))
  if (!is.na(proj.yr))
    p <- p + geom_vline(xintercept = proj.yr, linetype = "longdash")
  print(p)
  
}

