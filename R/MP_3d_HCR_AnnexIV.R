#-------------------------------------------------------------------------------
#                               annexIVHCR
#   HCR that implements the HCR described in Annex IV of (EC COM 2011 241)
#       rules 6 and 9.
#
#   Bnow = (Iy-1 + Iy-2)/2
#   Bref = (Iy-3 + Iy-4 + Iy-5)/3
#
#   TAC[y+1] = gamma*TAC[y]
#
#                   | 1 + beta      Bnow/Bref > 1 + alpha
#           gamma = | f(Bnow/Bref)  1-alpha < Bnow/Bref < 1 + alpha
#                   | 1 - beta      Bnow/Bref < 1 - alpha
#
#   Base Case: alpha = 0.2 and beta = 0.15
#
#   HCR 2 => f(Bnow/Bref)
#   HCR 4 => f(Bnow/Bref) = beta/alpha * (Bnow/Bref - 1) + 1
#
# * A single biomass index per stock.
#
# 13/09/2011 09:47:24
#-------------------------------------------------------------------------------

annexIVHCR <- function(indices, advice, advice.ctrl, year, stknm,...){

    Id <- indices[[stknm]][[1]]@index
    
    # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
    year.or <- year
    yrnm    <- dimnames(advice$TAC)[[2]][year]
    year    <- which(yrnm == dimnames(Id)[[2]])
    
    
    Bnow <- (yearSums(Id[,(year-1):(year-2)])/2)[drop=T] # [it]
    Bref <- (yearSums(Id[,(year-3):(year-5)])/3)[drop=T] # [it]
    Brat <- Bnow/Bref  [drop=T] # [it]
    
    alpha <- advice.ctrl[[stknm]][['ref.pts']]['alpha',]
    beta  <- advice.ctrl[[stknm]][['ref.pts']]['beta',]
    type  <- advice.ctrl[[stknm]][['type']]
    
    if(Brat > (1+alpha))
        gamma <- 1 + beta
    else{
        if(Brat < 1-alpha)
            gamma <- 1 - beta
        else{
            if(type == 2) gamma <- 1
            else{ # type == 4
                gamma <- beta/alpha * (Brat - 1) + 1
            }
        }
    }
    
    advice$TAC[stknm,year.or+1,] <- advice$TAC[stknm,year.or,]*gamma
    
    return(advice)
 }

