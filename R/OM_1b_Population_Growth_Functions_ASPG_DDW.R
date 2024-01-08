# Age structured population growth with densodependence
ASPG_DDW <- function(biols, SRs, fleets, year, season, stknm, biols.ctrl,...){
  
  # update wt of stknm
  ddw.model <- biols.ctrl[[stkm]][['DDW.model']]
  lfd       <- biols.ctrl[[stkm]][['LFD']]
  a         <- biols.ctrl[[stkm]][['a.lw']]
  lbins     <- as.numeric(colNames(lfd))
  
  B <- quantSums((biols[[stkm]]@wt*biols[[stkm]]@n)[,year-1])[drop=T]
    
  condF <- predict(LW_lm, data.frame(biomass = B))

  wy <- a*(lbins)^condF
  
  wt <- rowSums(sweep(lfd, 2, wy, "*"))
  
  biols[[stkm]]@wt[,year] <- wt
  
  # use normal ASPG to project the population
  
  res <- ASPG(biols, SRs, fleets, year, season, stknm, biols.ctrl,...)
  
  return(res)
} 


