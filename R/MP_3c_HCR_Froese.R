#-------------------------------------------------------------------------------
# HCR proposed in Froese et. al 2010 - Fish & Fisheries
# Biomass Based HCR.
#  Reference Points: Btrigger, Blim and Btarget.
#   The proposal in the paper is:
#           - Btrigger = 1*Bmsy
#           - Btarget  = 1.3*Bmsy
#           - Blim     = 0.5*Bmsy
#
#  - alpha_0 = 1, alpha_1 = 1.3 and alpha_2 = 0.5 are optional in this
#       implementation of the HCR.
#  - TAC advice depending on B in relation to BRP is:
#           - 0.91*MSY
#           - 0.91*f[beta]*MSY
#           - 0.
#   - beta = 0.91 is optional in this implementation of the HCR.
#
# 06/09/2011 12:11:51
#-------------------------------------------------------------------------------


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
    mult <- 1/(1-ref.pts['alpha_0',])*(-ref.pts['alpha_0',] + b.datyr/Brefs[2,])
    TAC  <- ifelse(b.pos == 0, 0, ifelse(b.pos == 1, ref.pts['MSY',]*ref.pts['beta',]*mult,ref.pts['MSY',]*ref.pts['beta',]))
    
    advice[['TAC']][stknm,year+1,,,,] <- TAC

    return(advice)
}

