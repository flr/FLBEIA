#-------------------------------------------------------------------------------
# HCR proposed in Froese et. al 2011 - Fish & Fisheries
# Biomass Based HCR.
#  Reference Points: Btrigger, Blim and Btarget.
#   The proposal in the paper is:
#           - Blim     = alpha_0 * Bmsy
#           - Btrigger = alpha_1 * Bmsy
#           - Btarget  = alpha_2 * Bmsy (not required for the rule)
#
#  - alpha_0 = 0.5, alpha_1 = 1 are optional in this
#       implementation of the HCR.
#
#  - TAC advice depending on B in relation to BRP is:
#           - 0.0                 , if B < Btrigger
#           - beta * MSY * f[Bmsy], if Btrigger <= B < Btarget
#           - MSY                 , if B >= Btarget
#   - beta = 0.91 is optional in this implementation of the HCR.
#   - f[Bmsy] = (B/Bmsy - alpha_0) / (alpha_1 - alpha_0)
#
# 06/09/2011 12:11:51
# Changed: 05/03/2019 - Sonia Sanchez (some corrections)
#-------------------------------------------------------------------------------
#' @rdname annualTAC
#' @aliases FroeseHCR

FroeseHCR <- function(stocks, advice, advice.ctrl, year, stknm,...){

    stk       <- stocks[[stknm]]
    ageStruct <- ifelse(dim(stk@m)[1] > 1, TRUE, FALSE)

    ref.pts  <- advice.ctrl[[stknm]]$ref.pts # matrix[5,it]  rows = Bmsy, MSY, alpha_0, alpha_1, beta
    
    iter     <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)

    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

    #  Calcuate where we are in relation to reference biomasses.
    Brefs <- sweep(ref.pts[c('alpha_0', 'alpha_1'), ,drop=F],2, ref.pts['Bmsy', ,drop=F], "*")

    # Last SSB (Age structured) OR Biomass (Aggregated) estimate
    if(ageStruct)
        b.datyr <- ssb(stk)[,year-1,drop = TRUE] # [it]
    else
        b.datyr <- (stk@stock.n*stk@stock.wt)[,year-1,drop = TRUE] # [it]
        
    # Find where the SSB (Age structured) OR Biomass (Aggregated) is in relation to reference points.
    b.pos <- apply(matrix(1:iter,1,iter),2, function(i) findInterval(b.datyr[i], Brefs[,i]))  # [it]

    # Calculate the TAC.
    mult <- 1/(ref.pts['alpha_1',] - ref.pts['alpha_0',]) * 
      (b.datyr/ref.pts['Bmsy',] - ref.pts['alpha_0',])
    TAC  <- ifelse(b.pos == 0, 0, ifelse(b.pos == 1, ref.pts['MSY',]*ref.pts['beta',]*mult,
                                         ref.pts['MSY',]*ref.pts['beta',]))
    
    advice[['TAC']][stknm,year+1,,,,] <- TAC

    return(advice)
}

