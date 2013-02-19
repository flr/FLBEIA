#-------------------------------------------------------------------------------
# HCR E - proposed rule for BoBane long term management plan
# SSB Based HCR.
#
#  - TAC advice depending on SSB in relation to BRP is:
#       Units = thousand tons
#           - 0      , ssb <= 24
#           - 7      , 24 < ssb < 33
#           - hr*ssb , ssb>=33
#     TAC<=33
#
# 04/06/2012 15:55:35
#-------------------------------------------------------------------------------


#! COMPROBAR QUE VALOR DE year HAY QUE USAR (no habría que considerar también season?)
aneHCRE <- function(stocks, advice, advice.ctrl, year, stknm,...){

    stk       <- stocks[[stknm]]
    ageStruct <- ifelse(dim(stk@m)[1] > 1, TRUE, FALSE)

    iter     <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)

    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]

    #  Calcuate where we are in relation to reference biomasses.
    Brefs <- c(0,24,33)

    # Last SSB (Age structured) OR Biomass (Aggregated) estimate
    if(ageStruct)
        b.datyr <- apply( stk@stock.n * stk@stock.wt * stk@mat * exp (-stk@m.spwn * stk@m - stk@harvest * stk@harvest.spwn),c(2,6),sum)[,year,drop = TRUE] # [it]
    else
        b.datyr <- (stk@stock.n*stk@stock.wt)[,year,drop = TRUE] # [it]
        
    # Find where the SSB (Age structured) OR Biomass (Aggregated) is in relation to reference points.
    b.pos <- apply(matrix(1:iter,1,iter),2, function(i) findInterval(b.datyr[i], Brefs))  # [it]

    # Calculate the TAC.
    TAC  <- ifelse(b.pos == 0, 0, ifelse(b.pos == 1, 7, 0.3*b.datyr))
    
    ifelse( TAC >33, 33, TAC)
    advice[['TAC']][stknm,year,,,,] <- TAC

    return(advice)
}

