#-------------------------------------------------------------------------------
#                               annexIVHCR
#   HCR that implements the ICES HCR for Category 3 DLS stocks
#       (see ICES CM 2012/ACOM 68: Category 3 - Method 3.2)
#
#   x : number of years in the survey average for current status
#   z : total number of years in the survey to be considered
#   
#   TAC[y+1] = k * TAC[y]
#      Bnow = (I_{y-1} + I_{y-x})/x
#      Bref = (I_{y-x+1} + ... + I_{y-z})/(z-x)
#          | 1-ucaplow, Brat <= 1-ucaplow
#      k = | Bnow/Bref, 1-ucaplow < Brat < 1+ucapupp
#          | 1+ucapupp, Brat >= 1+ucapupp
#
#   Base Case: x = 2 and z = 5
#
# * A single biomass index per stock.
#
# Author : 2019-05-28 14:35:31 - ssanchez
# Changed: 
#-------------------------------------------------------------------------------
#' @rdname annualTAC
#' @aliases IcesCat3HCR

IcesCat3HCR <- function(indices, advice, advice.ctrl, year, stknm,...){

    Idnm <- advice.ctrl[[stknm]][['index']]  # either the name or the position of the index in FLIndices object.
    Id <- indices[[stknm]][[Idnm]]@index
    
    ni <- dim(Id)[6]
    
    # tac.lag: used when in-year advice
    tac.lag <- advice.ctrl[[stknm]][['tac.lag']]
    
    # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
    year.or <- year
    yrnm    <- dimnames(advice$TAC)[[2]][year]
    year    <- which(yrnm == dimnames(Id)[[2]])
    
    if (!is.null(advice.ctrl[[stknm]][['inyr.idx']])) 
      if (advice.ctrl[[stknm]][['inyr.idx']]==TRUE)
        year <- year + 1
    
    # Reference years
    if (is.null(advice.ctrl[[stknm]][['refnyrs']])) {
      refnyrs <- c(now=2, past=3)
    } else
      refnyrs <- advice.ctrl[[stknm]][['refnyrs']]
    
    x <- refnyrs["now"]
    z <- x + refnyrs["past"]
    
    # Reference vaules (indices & catches)
    
    Bnow <- (yearSums(Id[,(year-1):(year-x)])/x)[drop=T]       # [it]
    Bref <- (yearSums(Id[,(year-x-1):(year-z)])/(z-x))[drop=T] # [it]
    Brat <- Bnow/Bref[drop=T] # [it]
    
    Cref <- advice$TAC[stknm,year.or-tac.lag,] # catches default (previous TAC)
    
    if (!is.null(advice.ctrl[[stknm]]$cref.year)) 
      if( yrnm==ac(advice.ctrl[[stknm]]$cref.year))
        Cref <- advice.ctrl[[stknm]]$cref.value # specific reference TAC set for current year
    
    # Uncertainty cap
    # - NULL  : no uncertainty cap
    # - value : value of the cap in %1
    
    ucaplow <- ifelse( is.null(advice.ctrl[[stknm]]$unccap.low), Inf, advice.ctrl[[stknm]]$unccap.low)
    ucapupp <- ifelse( is.null(advice.ctrl[[stknm]]$unccap.upp), Inf, advice.ctrl[[stknm]]$unccap.upp)
    
    Brat <- ifelse( Brat > 1+ucapupp, 1+ucapupp, ifelse(Brat < 1-ucaplow, 1-ucaplow, Brat))
    
    # Initial TAC
    
    TAC <- Cref * Brat
    
    # Precautionary buffer
    
    if (!is.null(advice.ctrl[[stknm]]$prec.buffer)) {
      
      # year: must be the name of the years or "all"
      
      if (!is.character(advice.ctrl[[stknm]]$prec.buffer$year))
        stop( paste("advice.ctrl[['",stknm,"']]$prec.buffer$year must be character vector.",sep=""))
      
      pbuf.yr <- advice.ctrl[[stknm]]$prec.buffer$year
      
      if (pbuf.yr=="all" | yrnm %in% pbuf.yr) {
        pbuf <- advice.ctrl[[stknm]]$prec.buffer$value # must be 0, in no precautionary buffer
        TAC <- TAC * (1-pbuf)
      }
      
    }
    
    advice$TAC[stknm,year.or+1-tac.lag,] <- TAC
    
    return(advice)
 }

