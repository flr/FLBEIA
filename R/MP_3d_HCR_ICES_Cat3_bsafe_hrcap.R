#-------------------------------------------------------------------------------
#   Alternative approach for annexIVHCR with a biomass safeguard and harvest rate caps
#   (from WKDLSSLS2019: Nicola Walker)
#
#   IcesCat3HCR_bsafe_hrcap implements the ICES HCR for Category 3 DLS stocks
#       (see ICES CM 2012/ACOM 68: Category 3 - Method 3.2)
#   with a biomass safeguard and hr limits
#
#   x : number of years in the survey average for current status
#   z : total number of years in the survey to be considered
#   
#   TAC[y+1] = k * Irat * TAC[y]
#      Bnow = (I_{y-1} + I_{y-x})/x
#      Bref = (I_{y-x+1} + ... + I_{y-z})/(z-x)
#      Brat = Bnow/Bref
#          | 1-ucaplow, Brat <= 1-ucaplow
#      k = | Brat     , 1-ucaplow < Brat < 1+ucapupp
#          | 1+ucapupp, Brat >= 1+ucapupp
#      Irat = min(1,I_y/I_trigger), where Itrigger can be minum observed index in the historical part
#                                   or 1.4 * minimum observed index in the historical part
#                                   or any other alternative, such as quantile(Ihist,0.05)
# Where:
#             | hrmin   , TAC[y+1]/Bnow <= hrmin
#  TAC[y+1] = | TAC[y+1], hrmin < TAC[y+1]/Bnow < hrmax
#             | hrmax   , TAC[y+1]/Bnow >= hrmax
#
#   Base Case: x = 2 and z = 5
#
# * A single biomass index per stock.
#
# Author : 2019-05-28 14:35:31 - ssanchez
# Changed: 2020-09-04 08:00:00 - libaibarriaga (adapted from IcesCat3HCR)
#-------------------------------------------------------------------------------
#' @rdname annualTAC
#' @aliases IcesCat3HCR_bsafe_hrcap

IcesCat3HCR_bsafe_hrcap <- function(indices, advice, advice.ctrl, year, stknm,...){

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
    
    # Reference index: [it]
    if ( is.null(advice.ctrl[[stknm]]$Itrigger)) {
      Itrigger <- rep(1,ni)
    } else
      Itrigger <- advice.ctrl[[stknm]]$Itrigger
    if (length(Itrigger)==1) Itrigger <- rep(Itrigger,ni) # [it]
    
    # Reference values (indices & catches)
    
    Bnow <- (yearSums(Id[,(year-1):(year-x)])/x)[drop=T]       # [it]
    Bref <- (yearSums(Id[,(year-x-1):(year-z)])/(z-x))[drop=T] # [it]
    Brat <- Bnow/Bref[drop=T]                                  # [it]
    
    Irat <- pmin(1, Bnow/Itrigger)                             # [it]
    
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
    
    TAC <- Cref * Brat * Irat
    
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

    # harvest rates caps
    # - NULL  : no harvest rate cap
    # - value : min and max value of the cap
    
    if(is.null(advice.ctrl[[stknm]]$hrmin)){
      hrmin <- rep(-Inf, ni)  
    }else{
      hrmin <- advice.ctrl[[stknm]]$hrmin
    }

    if(is.null(advice.ctrl[[stknm]]$hrmax)){
      hrmax <- rep(Inf, ni)  
    }else{
      hrmax <- advice.ctrl[[stknm]]$hrmax
    }
    
    TAC <- Bnow * ifelse( TAC/Bnow > hrmax, hrmax, ifelse( TAC/Bnow < hrmin, hrmin, TAC/Bnow))
    
    advice$TAC[stknm,year.or+1-tac.lag,] <- TAC
    
    return(advice)
 }

