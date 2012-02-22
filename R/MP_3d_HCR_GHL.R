#-------------------------------------------------------------------------------
#                               Greenland Halibut
#       HCR that implements the HCR used for the management of 
#       Greenland Halibut in NAFO
#
#   3 abundance indices.
#
#   - Calculate a linear model for each index. 
#   - Index the independent variable and the year the explanatory, used
#     the last 5 years for the regression.
#   - Extract the slopes and calculate their mean.
#   - TAC = TAC + (lambda*slope)
#         * lambda: alpha_0 if slope < 0
#         * lambda: alpha_1 if slope > 0
#
# 31/01/2012 09:01:12
#-------------------------------------------------------------------------------

ghlHCR <- function(indices, advice, advice.ctrl, year, stknm,...){

    it <- dim(indices[[1]][[1]]@index)[6]
    
    slopes <- matrix(NA,3,it)
    
    yrs.lm <- (year-5):(year-1) 
    yrs.lm <- dimnames(advice$TAC)[[2]][yrs.lm]
    
    for(id in 1:3){
        Id <- quantSums(indices[[stknm]][[id]]@index*indices[[stknm]][[id]]@catch.wt)
        for(i in 1:it){
            Idi <- c(Id[,yrs.lm,,,,i])
            slopes[id,i] <- coef(lm(Idi~factor(yrs.lm)))[2]
        }
    }
    
    slp <- apply(slopes,2,mean)
    
    lambda <- ifelse(slp < 0, advice.ctrl[[stknm]]$ref.pts['alpha_0',], advice.ctrl[[stknm]]$ref.pts['alpha_1',])
    mult   <- lambda*slp
    beta   <- advice.ctrl[[stknm]]$ref.pts['beta',]
    mult   <- ifelse(mult > (1+ beta), (1+ beta), ifelse(mult < (1-beta), 1-beta, mult))
    advice$TAC[stknm,year+1,] <- c(advice$TAC[stknm,year,])*mult
    
    return(advice)
 }

