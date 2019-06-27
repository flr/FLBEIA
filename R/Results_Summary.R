#------------------------------------------------------------------------------#
#    Auxiliary functions to summarize the results
#
# F_flbeia(obj)
# SSB_flbeia(obj)
# B_flbeia(obj)
# R_flbeia(obj)
# C_flbeia(obj);  L_flbeia(obj); D_flbeia(obj)
# summary_flbeia(obj)
#  obj = FLBEIA output.
#  
#   
# Dorleta Garcia
# Created: 30/01/2011 20:50:27 
# Changed: 30/01/2011 20:50:32
#------------------------------------------------------------------------------#

#' @title Biological summary functions
#' 
#' @description These functions return the biomass (B), fishing mortality (F),  spawning stock biomass (SSB), recruitment (R), catches (C), landings (L) and discards (D) indicators. 
#' 
#' @param  obj The output of the FLBEIA function.
#' @param  years The  years for which the indicators are extracted.
#' 
#' @return B_flbeia, F_flbeia... return an array with three dimensions (stock, year and iter).
#' The summary_flbeia function returns an array with 4 dimensions (stock, year, iter, indicator) with the value of all the indicators. 
#'         
#'
#' @details 
#' 
#' \itemize{
#'       \item{B_flbeia}{ this function computes SSB.}
#'       \item{F_flbeia}{ this function computes fishing mortiality.}
#'       \item{SSB_flbeia}{ this function computes spawning stock biomass by species.}
#'       \item{R_flbeia}{ this function computes recruitment by stock. If the stock is defined by age this function the recruiment is computed. ;
#'                        If the stock is follows a biomass dynamics, this function gives the growth.}
#'       \item{C_flbeia}{ this function computes catches by fleets and stock.} 
#'       \item{L_flbeia}{ this function computes landings by fleets and stock.}
#'       \item{D_flbeia}{ this function computes the discards by fleets and stock.}
#'       \item{summary_flbeia}{ this function computes recruitment, SSB, fishing mortality,
#'                              total biomass for all stocks; 
#'                              and catch and landings for all fleets by year and stock.}
#'      }     


#------------------------------------------------------------------------------#
# F_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
# @export

#' @rdname summary_flbeia
#' @aliases F_flbeia
#dga 2019/05/30: make it seasonal.
F_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    ns     <- dim(obj$biols[[1]]@n)[4]
    yrnms  <- years
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    
    res <- array(dim = c(length(stknms), ny,ns, it), dimnames = list(stock = stknms,  year = yrnms, season = ssnms,iter = 1:it))
    
    for(stk in stknms){
        # harvest: * if age structured calculate it from 'n'.
        #          * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
        na <- dim(obj$biols[[stk]]@n)[1]
        
        if(na == 1){
            # Catch:
            catch <- apply(catchStock(obj$fleets, stk),c(2,4,6), sum)[,years,drop = TRUE] # [ny,ns,it]
            B     <- (obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,years,,,drop= TRUE] # [ny,ns, it] , 1st season biomass
            res[stk,,,] <- (catch/B)
        }
        else{ 
            fbar_age <- ac(obj$biols[[stk]]@range[c('minfbar')]:obj$biols[[stk]]@range[c('maxfbar')])
            
            Dnms <- list(age = fbar_age, year = yrnms, season = ssnms, iter = 1:it)
            aux  <- array(dim = c(length(fbar_age), ny,ns,it), dimnames = Dnms)           
            
            n.  <- array(unitSums(obj$biols[[stk]]@n)[fbar_age,years,,,drop=T], dim = c(length(fbar_age),ny,ns,it), dimnames = Dnms)
            m.  <- array(seasonSums(unitMeans(obj$biols[[stk]]@m))[fbar_age,years,drop=T], dim = c(length(fbar_age),ny,ns,it), dimnames = Dnms)
            c.  <- array(apply(catchStock(obj$fleets, stk),c(1:2,4,6), sum)[fbar_age,years, drop = TRUE], dim = c(length(fbar_age),ny,ns,it), dimnames = Dnms)
        
            fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
        
            for(ss in ssnms){
            for(y in yrnms){
                for(a in fbar_age){
                    for(i in 1:it){
                      if(is.na(n.[a,y,ss,i])){ 
                        aux[a,y,ss,i] <- NA
                      }
                      else{
                        if(n.[a,y,ss,i] == 0) aux[a,y,ss,i] <- 0
                        else{
                          xx <- try(uniroot(fobj, lower = 0, upper = 1e6, n = n.[a,y,ss,i], m=m.[a,y,ss,i], c = c.[a,y,ss,i])$root, silent = TRUE)
                          aux[a,y,,i] <- ifelse(class(xx) == 'try-error', NA, xx)
                        }
                      }     
                    }}}}
            
           res[stk,,,] <- apply(aux,2:4,mean) 
        }
    }
    return(res)
}


#------------------------------------------------------------------------------#
# SSB_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
#' @aliases SSB_flbeia
SSB_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)
    n <- FLCore::n
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    ns     <- dim(obj$biols[[1]]@n)[4]
    yrnms  <- years
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    

    stknms <- names(obj$biols)

    
    res <- array(dim = c(length(stknms), ny,ns, it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){ # SSB in any of the seasons
      
      res[stk,,,] <- apply(unitSums(n(obj$biols[[stk]])*wt(obj$biols[[stk]])*fec(obj$biols[[stk]])*mat(obj$biols[[stk]])*
                                     exp(-spwn(obj$biols[[stk]])*m(obj$biols[[stk]])))[,years,,], c(2,4,6), sum, na.rm=TRUE)[drop=T]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# B_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
#' @aliases B_flbeia
B_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
 
    ns     <- dim(obj$biols[[1]]@n)[4]
    yrnms  <- years
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    

    stknms <- names(obj$biols)

    
    res <- array(dim = c(length(stknms), ny,ns,it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){ # B 1st season
        res[stk,,,] <- apply(unitSums(obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,years,,], c(2,4,6), sum,  na.rm=TRUE)[drop=T]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# R_flbeia(obj) :: res[stocks, years, it] 
# If age struc => recruitment.
# If biodyn    => growth.
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
#' @aliases R_flbeia
R_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)

    ns     <- dim(obj$biols[[1]]@n)[4]
    yrnms  <- years
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    
    

    stknms <- names(obj$biols)

    res <- array(0,dim = c(length(stknms), ny,ns,it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){
      if(dim(obj$biols[[stk]]@n)[[1]] > 1) for(ss in 1:ns) for(ss0 in 1:ss) res[stk,,ss,] <- res[stk,,ss,] + obj$biols[[stk]]@n[1,yrnms,ss0,ss0,drop=T]
      else{
            catch <- array(apply(catchStock(obj$fleets, stk),c(2,6), sum)[,years,drop = TRUE], dim = c(ny,ns,it))
            B     <- array((obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,years,,,drop= TRUE], dim = c(ny,ns,it))
            res[stk,-ny,,] <- B[-1,,] - B[-ny,,] + catch[-ny,,]
        }
    }
    return(res)
}


#------------------------------------------------------------------------------#
# C_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
#' @aliases C_flbeia
C_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)

    ns     <- dim(obj$biols[[1]]@n)[4]
    yrnms  <- years
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    
    stknms <- names(obj$biols)
    
    res <- array(dim = c(length(stknms), ny,ns,it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){ # B 1st season
        res[stk,,,] <- apply(catchWStock(obj$fleets, stk),c(2,4,6), sum)[,years,drop = TRUE] # [ny,it]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# L_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
#' @aliases L_flbeia
L_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    ns     <- dim(obj$biols[[1]]@n)[4]
    yrnms  <- years
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    

    stknms <- names(obj$biols)

    res <- array(dim = c(length(stknms), ny,ns,it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){ # B 1st season
        res[stk,,,] <- apply(landWStock(obj$fleets, stk),c(2,4,6), sum)[,years,drop = TRUE] # [ny,it]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# D_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
#' @aliases D_flbeia
D_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years
    ns     <- dim(obj$biols[[1]]@n)[4]
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    

    stknms <- names(obj$biols)

    res <- array(dim = c(length(stknms), ny,ns,it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){ # B 1st season
        res[stk,,,] <- apply(discWStock(obj$fleets, stk),c(2,4,6), sum)[,years,drop = TRUE] # [ny,it]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# SUMMARY OF FLBEIA output :: res[stocks, years, it, indicators] 
#------------------------------------------------------------------------------#

#' @rdname summary_flbeia
summary_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){

    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years) # dim(obj$biols[[1]]@n)[2]
    yrnms  <- years # dimnames(obj$biols[[1]]@n)[[2]]
    ns     <- dim(obj$biols[[1]]@n)[4]
    ssnms <- dimnames(obj$biols[[1]]@n)[[4]]
    
    res <- array(dim = c(length(stknms), ny,ns,it, 7), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it, 
                                                      indicators = c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards')))
    
    res[,,,,1] <- R_flbeia(obj,years)
    res[,,,,2] <- SSB_flbeia(obj,years)
    res[,,,,3] <- F_flbeia(obj,years)
    res[,,,,4] <- B_flbeia(obj,years)
    res[,,,,5] <- C_flbeia(obj,years)
    res[,,,,6] <- L_flbeia(obj,years)
    res[,,,,7] <- D_flbeia(obj,years)

    return(res)
    
}

#------------------------------------------------------------------------------#
# summary_flbeia(obj) :: res[stocks, years, it, indicators] 
#------------------------------------------------------------------------------#
#' Summary of the FLBEIA output 
#' 
#' Summarize the results of the simulation in data frames.
#'
#' @details 
#' 
#'\itemize{
#'      \item{advSum, advSumQ:} Data frame with the indicators related with the management advice (TAC). The indicators are:
#'              "catch", "discards", "discRat",  "landings", "quotaUpt" and "tac".              
#'      \item{bioSum, bioSumQ:} Data frame with the biological indicators. The indicators are: 
#'              "biomass", "catch", "catch.iyv", "discards",  "disc.iyv",  "f", "landings",  "land.iyv",  "rec" and      "ssb".
#'      \item{fltSum, fltSumQ:} Data frame with the indicators at fleet level. The indicators are:
#'              "capacity", "catch", "costs", "discards", "discRat", "effort",       
#'              "fcosts", "gva", "grossValue", "landings", "fep", "nVessels", "price", "grossSurplus",
#'              "quotaUpt", "salaries", "vcosts" and "profitability".
#'      \item{fltStkSum, fltStkSumQ:} Data frame with the indicators at fleet and stock level. The indicators are:
#'              "landings", "discards", "catch", "price",  "quotaUpt", "tacshare", "discRat" and  "quota".   
#'      \item{npv:} A data frame with the net present value per fleet over the selected range of years.
#'      \item{mtSum, mtSumQ:} Data frame with the indicators at metier level. The indicators are:
#'              "effshare", "effort", "grossValue" and "vcost".   
#'      \item{mtStkSum, mtStkSumQ:} Data frame with the indicators at fleet and metier level. The indicators are:
#'              "catch",  "discards", "discRat", "landings" and "price".
#'      \item{riskSum:} A data frame with the risk indicators. The indicators are:
#'              "pBlim", "pBpa" and "pPrflim".
#'      \item{vesselSum, vesselSumQ:} Data frame with the indicators at vessel level. The indicators are:
#'               "catch", "costs", "discards", "discRat", "effort",       
#'              "fcosts", "gva", "grossValue", "landings", "fep",  "price", "grossSurplus",
#'              "quotaUpt", "salaries", "vcosts" and "profitability".
#'      \item{vesselStkSum, vesselStkSumQ:} Data frame with the indicators at vessel and stock level. The indicators are:
#'              "landings", "discards", "catch", "price",  "quotaUpt", "tacshare", "discRat" and  "quota".   
#'      \item{summary_flbeia:} An array with four dimensions: stock, year, iteration, 
#'      indicator. The indicators are: recruitment, ssb, f, biomass, catch, landings and discards.
#'      \item{ecoSum_damara:} ecoSum built in the framework of Damara project.
#'}
#'      
#' The data frames       
#'
#'
#' @return A data frame with columns for scenario, year, stock, iter, indicator, value,...
#' 
#'  The data frames can be of wide or long format. In long format all the indicators are in the same column.
#'  There is one column, indicator, for the name of the indicator and a second one value for the numeric value of the indicator.
#'  In the wide format each of the indicators correspond with one column in the data frame. 
#'  The long format it is recommendable to work with ggplot2 functions for example while the wide format 
#'  it is more efficient for memory allocation and speed of computations. 
#'  
#'  The quantile version of the summaries, fooQ, returns the quantiles of the indicators. 
#'  In the long format as many columns as elements in prob are created. The name of the columns are
#'  the elements in prob preceded by a q. In the wide format for each of the indicators as many 
#'  columns as elements in prob are created. The names of the colums are the elements in prob preceded by
#'  q_name_of_the_indicator. 
#'
#' @inheritParams FLBEIA
#' @inheritParams summary_flbeia
#' @param flnms Names of the fleet for which the indicators will be calculated.
#' @param stknms Names of the stock for which the indicators will be calculated.
#' @param years the names of the years for which the indicators will be calculated. 
#' @param long logical. The data frame should be constructed using long or wide format? Default TRUE.
#' @param byyear logical. The indicators should be provided at season or year level? Default TRUE.
#' @param ssb_season If byyear = TRUE, the season in which ssb will be taken.
#' @param prob a numeric vector with the probabilities used to calculate the quantiles. 
#' @param scenario a character string with the name of the scenario corresponding with obj. Default bc.
#' @param Bpa named numeric vector with one element per stock in stknms. The precautionary approach stock spawning biomass used in riskSum function to calculate biological risk yearly.
#' @param Blim named numeric vector with one element per stock in stknms. The limit stock spawning biomass used in riskSum function to calculate biological risk yearly.
#' @param Prflim named numeric vector with one element per fleet in flnms. The limit profit level used in riskSum function to calculate economic risk yearly.
#' @param discF Discount rate.
#' @param y0 character. Reference year.

#' @examples
#'\dontrun{
#'
#' library(FLBEIA)
#'
#' # Apply the summary functions to the examples runs in FLBEIA help page.
#' # Test the different arguments in summary function.
#' 
#' data(res_flbeia)
# 
#' #------------------------------------------------
#' # Example One: One stock, one fleet, one iter.
#' #------------------------------------------------
#' oneRes_bio    <- bioSum(oneRes)
#' oneRes$fleets[[1]] <- setUnitsNA(oneRes$fleets[[1]]) 
#' oneRes_flt    <- fltSum(oneRes)
#' oneRes_fltStk <- fltStkSum(oneRes)
#' oneRes_mt     <- mtSum(oneRes)
#' oneRes_mtStk  <- mtStkSum(oneRes)
#' oneRes_adv    <- advSum(oneRes)
#' 
#' head(oneRes_bio)
#' head(oneRes_flt)
#' head(oneRes_fltStk)
#' head(oneRes_mt)
#' head(oneRes_mtStk)
#' head(oneRes_adv)
#' 
#' oneRes_bioQ    <- bioSumQ(oneRes_bio)
#' oneRes_fltQ    <- fltSumQ(oneRes_flt)
#' oneRes_fltStkQ <- fltStkSumQ(oneRes_fltStk)
#' oneRes_mtQ     <- mtSumQ(oneRes_mt)
#' oneRes_mtStkQ  <- mtStkSumQ(oneRes_mtStk)
#' oneRes_advQ    <- advSumQ(oneRes_adv)
#' 
#' head(oneRes_bioQ)
#' head(oneRes_fltQ)
#' head(oneRes_fltStkQ)
#' head(oneRes_mtQ)
#' head(oneRes_mtStkQ)
#' head(oneRes_advQ)
#' 
#' # Wide format
#' oneRes_bio    <- bioSum(oneRes, long = TRUE, years = ac(2016:2020))
#' oneRes_flt    <- fltSum(oneRes, long = TRUE, years = ac(2016:2020))
#' oneRes_fltStk <- fltStkSum(oneRes, long = TRUE, years = ac(2016:2020))
#' oneRes_mt     <- mtSum(oneRes, long = TRUE, years = ac(2016:2020))
#' oneRes_mtStk  <- mtStkSum(oneRes, long = TRUE, years = ac(2016:2020))
#' oneRes_adv    <- advSum(oneRes, long = TRUE, years = ac(2016:2020))
#' 
#' head(oneRes_bio)
#' head(oneRes_flt)
#' head(oneRes_fltStk)
#' head(oneRes_mt)
#' head(oneRes_mtStk)
#' head(oneRes_adv)
#' 
#' oneRes_bioQ    <- bioSumQ(oneRes_bio)
#' oneRes_fltQ    <- fltSumQ(oneRes_flt)
#' oneRes_fltStkQ <- fltStkSumQ(oneRes_fltStk)
#' oneRes_mtQ     <- mtSumQ(oneRes_mt)
#' oneRes_mtStkQ  <- mtStkSumQ(oneRes_mtStk)
#' oneRes_advQ    <- advSumQ(oneRes_adv)
#' 
#' head(oneRes_bio)
#' head(oneRes_flt)
#' head(oneRes_fltStk)
#' head(oneRes_mt)
#' head(oneRes_mtStk)
#' head(oneRes_adv)
#' 
#' # Wide format with seasonal disaggregation. No seasonal disagregation available for
#' #  adv summaries.
#' 
#' oneRes_bio    <- bioSum(oneRes, long = FALSE, byyear = FALSE) # Biol summary is only by year.
#' oneRes_flt    <- fltSum(oneRes, long = FALSE, byyear = FALSE)
#' oneRes_fltStk <- fltStkSum(oneRes, long = FALSE, byyear = FALSE)
#' oneRes_mt     <- mtSum(oneRes, long = FALSE, byyear = FALSE)
#' oneRes_mtStk  <- mtStkSum(oneRes, long = FALSE, byyear = FALSE)
#' oneRes_adv    <- advSum(oneRes, long = FALSE) # Advice summary is only by year.
#' 
#' oneRes_bioQ    <- bioSumQ(oneRes_bio)
#' oneRes_fltQ    <- fltSumQ(oneRes_flt)
#' oneRes_fltStkQ <- fltStkSumQ(oneRes_fltStk)
#' oneRes_mtQ     <- mtSumQ(oneRes_mt)
#' oneRes_mtStkQ  <- mtStkSumQ(oneRes_mtStk)
#' oneRes_advQ    <- advSumQ(oneRes_adv)
#' # #  # Long format and seasonaloneRes_bio    <- bioSum(oneRes, long = TRUE) # Biol summary is only by year.
#' oneRes_flt    <- fltSum(oneRes, long = TRUE, byyear = FALSE)
#' oneRes_fltStk <- fltStkSum(oneRes, long = TRUE, byyear = FALSE)
#' oneRes_mt     <- mtSum(oneRes, long = TRUE, byyear = FALSE)
#' oneRes_mtStk  <- mtStkSum(oneRes, long = TRUE, byyear = FALSE)
#' oneRes_adv    <- advSum(oneRes, long = TRUE) # Advice summary is only by year.
#' 
#' oneRes_bioQ    <- bioSumQ(oneRes_bio)
#' oneRes_fltQ    <- fltSumQ(oneRes_flt)
#' oneRes_fltStkQ <- fltStkSumQ(oneRes_fltStk)
#' oneRes_mtQ     <- mtSumQ(oneRes_mt)
#' oneRes_mtStkQ  <- mtStkSumQ(oneRes_mtStk)
#' oneRes_advQ    <- advSumQ(oneRes_adv)
#' 
#' 
#' #------------------------------------------------
#' # Example OneIt: As one but with iterations.
#' #------------------------------------------------
#' oneItRes_bio    <- bioSum(oneItRes, scenario = 'with_iters')
#' oneItRes$fleets[[1]] <- setUnitsNA(oneItRes$fleets[[1]])
#' oneItRes_flt    <- fltSum(oneItRes, scenario = 'with_iters')
#' oneItRes_fltStk <- fltStkSum(oneItRes, scenario = 'with_iters')
#' oneItRes_mt     <- mtSum(oneItRes, scenario = 'with_iters')
#' oneItRes_mtStk  <- mtStkSum(oneItRes, scenario = 'with_iters')
#' oneItRes_adv    <- advSum(oneItRes, scenario = 'with_iters')
#' 
#' oneItRes_bioQ    <- bioSumQ(oneItRes_bio)
#' oneItRes_fltQ    <- fltSumQ(oneItRes_flt)
#' oneItRes_fltStkQ <- fltStkSumQ(oneItRes_fltStk)
#' oneItRes_mtQ     <- mtSumQ(oneItRes_mt)
#' oneItRes_mtStkQ  <- mtStkSumQ(oneItRes_mtStk)
#' oneItRes_advQ    <- advSumQ(oneItRes_adv)
#' 
#' oneItRes_bio    <- bioSum(oneItRes, long = FALSE, years = ac(2016:2020))
#' oneItRes_flt    <- fltSum(oneItRes, long = FALSE, years = ac(2016:2020))
#' oneItRes_fltStk <- fltStkSum(oneItRes, long = FALSE, years = ac(2016:2020))
#' oneItRes_mt     <- mtSum(oneItRes, long = FALSE, years = ac(2016:2020))
#' oneItRes_mtStk  <- mtStkSum(oneItRes, long = FALSE, years = ac(2016:2020))
#' oneItRes_adv    <- advSum(oneItRes, long = FALSE, years = ac(2016:2020))
#' 
#'
#' oneItRes_bioQ    <- bioSumQ(oneItRes_bio)
#' oneItRes_fltQ    <- fltSumQ(oneItRes_flt)
#' oneItRes_fltStkQ <- fltStkSumQ(oneItRes_fltStk)
#' oneItRes_mtQ     <- mtSumQ(oneItRes_mt)
#' oneItRes_mtStkQ  <- mtStkSumQ(oneItRes_mtStk)
#' oneItRes_advQ    <- advSumQ(oneItRes_adv)
#' 
#' 
#' oneItRes_bio    <- bioSum(oneItRes, long = FALSE) # Biol summary is only by year.
#' oneItRes_flt    <- fltSum(oneItRes, long = FALSE, byyear = FALSE)
#' oneItRes_fltStk <- fltStkSum(oneItRes, long = FALSE, byyear = FALSE)
#' oneItRes_mt     <- mtSum(oneItRes, long = FALSE, byyear = FALSE)
#' oneItRes_mtStk  <- mtStkSum(oneItRes, long = FALSE, byyear = FALSE)
#' oneItRes_adv    <- advSum(oneItRes, long = FALSE) # Advice summary is only by year.
#' 
#' oneItRes_bioQ    <- bioSumQ(oneItRes_bio)
#' oneItRes_fltQ    <- fltSumQ(oneItRes_flt)
#' oneItRes_fltStkQ <- fltStkSumQ(oneItRes_fltStk)
#' oneItRes_mtQ     <- mtSumQ(oneItRes_mt)
#' oneItRes_mtStkQ  <- mtStkSumQ(oneItRes_mtStk)
#' oneItRes_advQ    <- advSumQ(oneItRes_adv)
#' 
#' 
#' oneItRes_bio    <- bioSum(oneItRes, long = TRUE) # Biol summary is only by year.
#' oneItRes_flt    <- fltSum(oneItRes, long = TRUE, byyear = FALSE)
#' oneItRes_fltStk <- fltStkSum(oneItRes, long = TRUE, byyear = FALSE)
#' oneItRes_mt     <- mtSum(oneItRes, long = TRUE, byyear = FALSE)
#' oneItRes_mtStk  <- mtStkSum(oneItRes, long = TRUE, byyear = FALSE)
#' oneItRes_adv    <- advSum(oneItRes, long = TRUE) # Advice summary is only by year.
#' 
#' oneItRes_bioQ    <- bioSumQ(oneItRes_bio)
#' oneItRes_fltQ    <- fltSumQ(oneItRes_flt)
#' oneItRes_fltStkQ <- fltStkSumQ(oneItRes_fltStk)
#' oneItRes_mtQ     <- mtSumQ(oneItRes_mt)
#' oneItRes_mtStkQ  <- mtStkSumQ(oneItRes_mtStk)
#' oneItRes_advQ    <- advSumQ(oneItRes_adv)
#' 
#' oneItRes_risk <- riskSum( oneItRes, Bpa = c(stk1= 900), Blim = c(stk1 = 600),
#'                           Prflim = c(fl1 = 0), scenario = 'alternative')
#' 
#' oneItRes_npv  <- npv(oneItRes, y0 = '2014')
#' 
#' #------------------------------------------------
#' # Example Multi: Two stock, two fleet, four iters.
#' #------------------------------------------------
#' multiRes_bio    <- bioSum(multiRes)
#' multiRes$fleets <- FLFleetsExt(lapply(multiRes$fleets, function(x) setUnitsNA(x)))
#' multiRes_flt    <- fltSum(multiRes)
#' multiRes_fltStk <- fltStkSum(multiRes)
#' multiRes_mt     <- mtSum(multiRes)
#' multiRes_mtStk  <- mtStkSum(multiRes)
#' multiRes_adv    <- advSum(multiRes)
#' 
#' multiRes_bioQ    <- bioSumQ(multiRes_bio)
#' multiRes_fltQ    <- fltSumQ(multiRes_flt)
#' multiRes_fltStkQ <- fltStkSumQ(multiRes_fltStk)
#' multiRes_mtQ     <- mtSumQ(multiRes_mt)
#' multiRes_mtStkQ  <- mtStkSumQ(multiRes_mtStk)
#' multiRes_advQ    <- advSumQ(multiRes_adv)
#' 
#' multiRes_bio    <- bioSum(multiRes, long = FALSE, years = ac(2016:2020))
#' multiRes_flt    <- fltSum(multiRes, long = FALSE, years = ac(2016:2020))
#' multiRes_fltStk <- fltStkSum(multiRes, long = FALSE, years = ac(2016:2020))
#' multiRes_mt     <- mtSum(multiRes, long = FALSE, years = ac(2016:2020))
#' multiRes_mtStk  <- mtStkSum(multiRes, long = FALSE, years = ac(2016:2020))
#' multiRes_adv    <- advSum(multiRes, long = FALSE, years = ac(2016:2020))
#' 
#' 
#' multiRes_bioQ    <- bioSumQ(multiRes_bio)
#' multiRes_fltQ    <- fltSumQ(multiRes_flt)
#' multiRes_fltStkQ <- fltStkSumQ(multiRes_fltStk)
#' multiRes_mtQ     <- mtSumQ(multiRes_mt)
#' multiRes_mtStkQ  <- mtStkSumQ(multiRes_mtStk)
#' multiRes_advQ    <- advSumQ(multiRes_adv)
#' 
#' 
#' multiRes_bio    <- bioSum(multiRes, long = FALSE, byyear = FALSE)
#' multiRes_flt    <- fltSum(multiRes, long = FALSE, byyear = FALSE)
#' multiRes_fltStk <- fltStkSum(multiRes, long = FALSE, byyear = FALSE)
#' multiRes_mt     <- mtSum(multiRes, long = FALSE, byyear = FALSE)
#' multiRes_mtStk  <- mtStkSum(multiRes, long = FALSE, byyear = FALSE)
#' multiRes_adv    <- advSum(multiRes, long = FALSE) # Advice summary is only by year.
#' 
#' multiRes_bioQ    <- bioSumQ(multiRes_bio)
#' multiRes_fltQ    <- fltSumQ(multiRes_flt)
#' multiRes_fltStkQ <- fltStkSumQ(multiRes_fltStk)
#' multiRes_mtQ     <- mtSumQ(multiRes_mt)
#' multiRes_mtStkQ  <- mtStkSumQ(multiRes_mtStk)
#' multiRes_advQ    <- advSumQ(multiRes_adv)
#' 
#' 
#' multiRes_bio    <- bioSum(multiRes, long = TRUE, byyear = FALSE)
#' multiRes_flt    <- fltSum(multiRes, long = TRUE, byyear = FALSE)
#' multiRes_fltStk <- fltStkSum(multiRes, long = TRUE, byyear = FALSE)
#' multiRes_mt     <- mtSum(multiRes, long = TRUE, byyear = FALSE)
#' multiRes_mtStk  <- mtStkSum(multiRes, long = TRUE, byyear = FALSE)
#' multiRes_adv    <- advSum(multiRes, long = TRUE) # Advice summary is only by year.
#' 
#' multiRes_bioQ    <- bioSumQ(multiRes_bio)
#' multiRes_fltQ    <- fltSumQ(multiRes_flt)
#' multiRes_fltStkQ <- fltStkSumQ(multiRes_fltStk)
#' multiRes_mtQ     <- mtSumQ(multiRes_mt)
#' multiRes_mtStkQ  <- mtStkSumQ(multiRes_mtStk)
#' multiRes_advQ    <- advSumQ(multiRes_adv)
#' 
#' multiRes_npv  <- npv(multiRes, y0 = '2014')
#' risk_multiRes <- riskSum( multiRes, Bpa = c(stk1= 135000, stk2 = 124000),
#'                           Blim = c(stk1= 96000, stk2 = 89000), Prflim = c(fl1 = 0, fl2 = 0),
#'                           scenario = 'alternative')
#' 
#' }


#------------------------------------------------------------------------------#
# bioSum :: data.frame[scenario, year, stock, iter, ||,||
#        rec, ssb, f, biomass, catch, landings, discards, land.iyv, disc.iyv, catch.iyv]
#------------------------------------------------------------------------------#
bioSum <- function(obj, stknms = 'all', years = dimnames(obj$biols[[1]]@n)$year, long = FALSE, scenario = 'bc', byyear = TRUE, ssb_season = NULL){
 
  if(stknms == 'all') stknms <- names(obj$biols)  
  
  if(byyear == TRUE & is.null(ssb_season)) ssb_season <- dimnames(obj$biols[[1]]@n)[[4]][1]
  
  xx <- summary_flbeia(obj) # array: stk x year x season x iter x indicator
  dat <- array2df(xx, label.x="value")
  
  # Wide format 
  res <- dat %>% mutate(scenario = scenario) %>% tidyr::spread(key=indicators, value=value) %>% 
    filter(year %in% years & stock %in% stknms) %>% arrange(year) %>%
    dplyr::group_by(stock, year, season, iter, scenario) %>% mutate(catch.iyv = catch/lag(catch), land.iyv = landings/lag(landings), 
                                                   disc.iyv = discards/lag(discards)) 
  
  # year or seasonal?
  if(byyear == TRUE){
    if(length(unique(res$season)) >1){
      # indicators that are summ up over the seasons
      res1 <- res %>% dplyr::group_by(stock, year, iter, scenario)  %>% summarise_at(c('rec','f','catch', 'landings', 'discards'),'sum')
      # ssb user selects the season
      res2 <- res %>%  filter(season == ssb_season) %>% ungroup() %>% select(stock, year, iter, scenario, ssb)
      # biomass the first season
      res3 <- res %>%  filter(season == dimnames(obj$biols[[1]]@n)[[4]][1]) %>% ungroup() %>% select(stock, year, iter, scenario, biomass)
      
      res <- full_join(res1, res2, by = c('stock', 'year', 'iter', 'scenario'))
      res <- full_join(res, res3, by = c('stock', 'year', 'iter', 'scenario'))
      res <- res %>% arrange(year) %>%  dplyr::group_by(stock, year, iter, scenario) %>%  mutate(catch.iyv = catch/lag(catch), land.iyv = landings/lag(landings),
                                                         disc.iyv = discards/lag(discards)) 
    }
    else{
      res <- res %>% ungroup() %>% select(-season)
    }
  }
  
  # reshaping this to the long format
  if(long == TRUE)
  res <- res %>% gather(key=indicator, value=value, biomass, catch, discards, f, landings, rec, ssb, catch.iyv, land.iyv, disc.iyv)
  
  return(res)
}

#' @rdname bioSum
#' @aliases bioSumQ
bioSumQ <- function(obj,  probs = c(0.95,0.5,0.05)){

  p_names <- paste("q",ifelse(nchar(substr(probs,3, nchar(probs)))==1, 
                              paste(substr(probs,3, nchar(probs)), 0, sep = ""), 
                              substr(probs,3, nchar(probs))), sep = "")
  
  if(dim(obj)[2] <= 7){ # the object is in long format
    
    if('season' %in% names(obj))
      res <- obj %>% dplyr::group_by(stock, year, season, scenario, indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(value, probs=probs, na.rm = TRUE))) %>% 
        unnest %>% tidyr::spread(key=quantiles, value=value)
    
    else
      
      res <- obj %>% dplyr::group_by(stock, year, scenario, indicator) %>% 
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(value, probs=probs, na.rm = TRUE))) %>% 
        unnest %>% tidyr::spread(key=quantiles, value=value)
    
  }
  else{
    
    p_funs <- purrr::map(probs, ~purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    if('season' %in% names(obj)){
      
      res <- obj %>% dplyr::group_by(stock, year, season, scenario) %>%  
        summarise_at(c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards', 'catch.iyv', 'land.iyv', 'disc.iyv'),
                     .funs =  p_funs)
    }
    else{
      
      res <- obj %>% dplyr::group_by(stock, year, scenario) %>%  
        summarise_at(c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards', 'catch.iyv', 'land.iyv', 'disc.iyv'),
                     .funs =  p_funs)
      
    }
    
  }
  
  return(res)
  }


#------------------------------------------------------------------------------#
# fltSum :: data.frame[scenario, year, season, fleet, iter, ||,|| 
#        capacity, catch, costs, discards, discRat, effort, fcosts, gva, grossValue, 
#        landings, fep, nVessels, price, grossSurplus, quotaUpt, salaries, 
#        vcosts, profitability]
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases fltSum
#' @param InterestRate Capital oportunity cost rate.
fltSum <- function (obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = TRUE, InterestRate = 0.03,scenario = 'bc')
{
  fleets <- obj$fleets
  covars <- obj$covars
  
#  fleets <- lapply(fleets, setUnitsNA)
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')

    if (flnms[1] == "all") flnms <- names(fleets)
    
    Dim <- dim(fleets[[1]]@effort[,years , ])[c(2, 4, 6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years , ])
    n <- prod(Dim) * length(flnms)
    
    if(is.null(covars$Depreciation)){ 
      covars$Depreciation <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, , ])[2:6]))
      dimnames(covars$Depreciation) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, , ])[2:6])
    }
    if(is.null(covars$Salaries)){ 
      covars$Salaries <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, , ])[2:6]))
      dimnames(covars$Salaries) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, , ])[2:6])
    }
    if(is.null(covars$MaxDays)){ 
      covars$MaxDays <- FLQuant(365/dim(fleets[[1]]@effort[, years, ])[4], dim = c(length(fleets),dim(fleets[[1]]@effort[, , ])[2:6]))
      dimnames(covars$MaxDays) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, , ])[2:6])
    }
    if(is.null(covars$NumbVessels)){ 
      covars$NumbVessels <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, , ])[2:6]))
      dimnames(covars$NumbVessels) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, , ])[2:6])
    }
    if(is.null(covars$CapitalCost)){ 
      covars$CapitalCost <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, , ])[2:6]))
      dimnames(covars$CapitalCost) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, , ])[2:6])
    }
    
    
    if(byyear == F){
      res <- data.frame(year = rep(years, prod(Dim[2:3]) * length(flnms)),
        season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3] *
            length(flnms)), fleet = rep(flnms, each = prod(Dim)),
        iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(flnms)),
        capacity = numeric(n), 
        catch = numeric(n),
        costs = numeric(n), 
        discards = numeric(n),
        discRat = numeric(n),
        effort = numeric(n),
        fcosts = numeric(n),
 #       gcf = numeric(n),
        gva = numeric(n),
        grossValue = numeric(n), 
        landings = numeric(n),
        fep = numeric(n), 
        NetProfit = numeric(n),
        nVessels = numeric(n), 
        price = numeric(n), 
        grossSurplus = numeric(n),
        quotaUpt = numeric(n), 
        salaries = numeric(n), 
        vcosts   = numeric(n),
       stringsAsFactors = FALSE)
    
    k <- 1
    for (f in flnms) {
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(catchWStock.f(fl, x))))
        res[k:(k + prod(Dim) - 1), "catch"] <- c(Reduce('+',temp)[,years])
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(landWStock.f(fl, x))))
        res[k:(k + prod(Dim) - 1), "landings"] <- c(Reduce('+',temp)[,years])
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(discWStock.f(fl, x))))
        res[k:(k + prod(Dim) - 1), "discards"] <- c(Reduce('+',temp)[,years])
        
        res[k:(k + prod(Dim) - 1), "discRat"] <- res[k:(k + prod(Dim) - 1), "discards"]/res[k:(k + prod(Dim) - 1), "catch"]
          
        res[k:(k + prod(Dim) - 1), "capacity"] <- c(fl@capacity[,years, ])
        
        res[k:(k + prod(Dim) - 1), "effort"] <- c(fl@effort[,years, ])
        
        res[k:(k + prod(Dim) - 1), "fcosts"] <- c(totfcost_flbeia(fl, covars, f)[,years, ])
        
        res[k:(k + prod(Dim) - 1), "vcosts"] <- c(totvcost_flbeia(fl)[,years, ])
        
        res[k:(k + prod(Dim) - 1), "costs"] <- c(costs_flbeia(fl, covars, f)[,years, ])
       
        res[k:(k + prod(Dim) - 1), "grossValue"] <- c(revenue_flbeia(fl)[,years, ]) 
        
        res[k:(k + prod(Dim) - 1), "grossSurplus"] <- c(revenue_flbeia(fl)[,years, ]) - res[k:(k + prod(Dim) - 1), "costs"]
        
        res[k:(k + prod(Dim) - 1), "price"] <- res[k:(k + prod(Dim) - 1), "grossValue"] / res[k:(k + prod(Dim) - 1), "landings"]
         
        res[k:(k + prod(Dim) - 1), "salaries"] <- c(fl@crewshare[,years,]*revenue_flbeia(fl)[,years, ] + covars[['Salaries']][f,years])
        
        res[k:(k + prod(Dim) - 1), "gva"] <- res[k:(k + prod(Dim) - 1), "grossValue"] -  res[k:(k + prod(Dim) - 1), "costs"] + res[k:(k + prod(Dim) - 1), "salaries"]
          
        
        res[k:(k + prod(Dim) - 1), "profitability"] <- res[k:(k + prod(Dim) - 1), "grossSurplus"]/res[k:(k + prod(Dim) - 1), "grossValue"]
 
        res[k:(k + prod(Dim) - 1), "nVessels"]  <- c(covars[['NumbVessels']][f,years])
          
        res[k:(k + prod(Dim) - 1), "fep"] <-  res[k:(k + prod(Dim) - 1), "grossSurplus"] - c(covars[['Depreciation']][f,years]*covars[['NumbVessels']][f,years])
        
        res[k:(k + prod(Dim) - 1), "netProfit"] <-  res[k:(k + prod(Dim) - 1), "fep"] - c(covars[['CapitalCost']][f,years]*InterestRate*covars[['NumbVessels']][f,years])
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(catchWStock.f(fl, x))))
        temp <- Reduce('+',temp)[,years]
        totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,years]*obj$advice$TAC[x,years]))
        if(dim(temp)[4] > 1) {res[k:(k + prod(Dim) - 1), "quotaUpt"] <- c(sweep(temp, c(1:3,5:6),totTAC/dim(temp)[4], "/"))}
        else{       res[k:(k + prod(Dim) - 1), "quotaUpt"] <- c(temp/totTAC)}
        
        k <- k + prod(Dim)
    }}
    else{
      n <- prod(Dim[-2]) * length(flnms)
      
      res <- data.frame(year = rep(years, prod(Dim[3]) * length(flnms)),
                        fleet = rep(flnms, each = prod(Dim[-2])),
                        iter = rep(rep(1:Dim[3], each = prod(Dim[1])), length(flnms)),
                        capacity = numeric(n), 
                        catch = numeric(n),
                        costs = numeric(n), 
                        discards = numeric(n),
                        discRat = numeric(n),
                        effort = numeric(n),
                        fcosts = numeric(n),
                        #       gcf = numeric(n),
                        gva = numeric(n),
                        grossValue = numeric(n), 
                        landings = numeric(n),
                        fep = numeric(n), 
                        nVessels = numeric(n), 
                        price = numeric(n), 
                        grossSurplus = numeric(n),
                        netProfit = numeric(n),
                        quotaUpt = numeric(n), 
                        salaries = numeric(n), 
                        vcosts   = numeric(n),
                        stringsAsFactors = FALSE)
      
      k <- 1
      for (f in flnms) {
   #     print(f)
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
        res[k:(k + prod(Dim[-2]) - 1), "catch"] <- c(Reduce('+',temp)[, years, ])
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(landWStock.f(fl, x)))))
        res[k:(k + prod(Dim[-2]) - 1), "landings"] <- c(Reduce('+',temp)[, years, ])
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(discWStock.f(fl, x)))))
        res[k:(k + prod(Dim[-2]) - 1), "discards"] <- c(Reduce('+',temp)[, years, ])
        
        res[k:(k + prod(Dim[-2]) - 1), "discRat"] <- res[k:(k + prod(Dim[-2]) - 1), "discards"]/res[k:(k + prod(Dim[-2]) - 1), "catch"]
        
        res[k:(k + prod(Dim[-2]) - 1), "capacity"] <- c(seasonSums(fl@capacity[,years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "effort"] <- c(seasonSums(fl@effort[,years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "fcosts"] <- c(seasonSums(totfcost_flbeia(fl, covars, f)[,years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "vcosts"] <- c(seasonSums(totvcost_flbeia(fl)[,years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "costs"] <- c(seasonSums(costs_flbeia(fl, covars, f)[,years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "grossValue"] <- c(seasonSums(revenue_flbeia(fl)[,years, ])) 
        
        res[k:(k + prod(Dim[-2]) - 1), "price"] <- res[k:(k + prod(Dim[-2]) - 1), "grossValue"] / res[k:(k + prod(Dim[-2]) - 1), "landings"]
        
        res[k:(k + prod(Dim[-2]) - 1), "grossSurplus"] <- c(seasonSums(revenue_flbeia(fl)[,years, ])) - res[k:(k + prod(Dim[-2]) - 1), "costs"]
        
        res[k:(k + prod(Dim[-2]) - 1), "salaries"] <- c(seasonSums(fl@crewshare[,years,]*revenue_flbeia(fl)[,years, ] + covars[['Salaries']][f,years]))
        
        res[k:(k + prod(Dim[-2]) - 1), "gva"] <- res[k:(k + prod(Dim[-2]) - 1), "grossValue"] - res[k:(k + prod(Dim[-2]) - 1), "costs"] + res[k:(k + prod(Dim[-2]) - 1), "salaries"]
        
        res[k:(k + prod(Dim[-2]) - 1), "profitability"] <- res[k:(k + prod(Dim[-2]) - 1), "grossSurplus"]/res[k:(k + prod(Dim[-2]) - 1), "grossValue"]
        
        res[k:(k + prod(Dim[-2]) - 1), "nVessels"]  <- c(seasonMeans(covars[['NumbVessels']][f, years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "fep"] <- c(seasonSums(revenue_flbeia(fl)[,years, ] -  costs_flbeia(fl, covars, f)[,years, ] - covars[['Depreciation']][f,years]*covars[['NumbVessels']][f,years]))
        
        res[k:(k + prod(Dim[-2]) - 1), "netProfit"] <- c(seasonSums(revenue_flbeia(fl)[,years, ] -  costs_flbeia(fl, covars, f)[,years, ] - covars[['Depreciation']][f,years]*covars[['NumbVessels']][f,years] - covars[['CapitalCost']][f,years]*InterestRate*covars[['NumbVessels']][f,years]))
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
        temp <- Reduce('+',temp)[, years, ]
        totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,years]*obj$advice$TAC[x,years]))
        res[k:(k + prod(Dim[-2]) - 1), "quotaUpt"] <- c(temp/totTAC)
        
        k <- k + prod(Dim[-2])
    }}
    
    if(long == TRUE){ # transform res into long format
      r1 <- ifelse(byyear == TRUE, 4,5)
      r2 <- ifelse(byyear == TRUE, 22,23)
      
      names(res)[r1:r2] <- paste('indicator',names(res)[r1:r2], sep = "_")
      res <- reshape(res, direction = 'long', varying = r1:r2, sep = "_")[,1:(r1+1)]
      rownames(res) <- 1:dim(res)[1]
      names(res)[!(names(res) %in% c("year", "season", "fleet", "iter", "time"))] <- 'value'
      names(res)[names(res) == 'time'] <- 'indicator'
      if('season' %in% names(res)) res <- res[,c('year','season', 'fleet', 'iter', 'indicator', 'value')]
      else res <- res[,c('year', 'fleet', 'iter', 'indicator', 'value')]
      # res[, r1:(r1+1)] <-  res[, (r1+1):r1]
      # names(res)[r1:(r1+1)] <- c( 'indicator', 'value') 
    }
  
   res <- cbind(scenario = scenario, res)
    return(res)
}
    
#' @rdname bioSum
#' @aliases fltSumQ
fltSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + indicator + year + scenario, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:4], data.frame(res[,5]))
    
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
      names(res)[5:(5+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + indicator + year + scenario + season, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[6:(6+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(capacity = obj$capacity,      catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                          discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                          grossValue  = obj$grossValue,         landings  = obj$landings,  fep  = obj$fep, netProfit = obj$netProfit, nVessels  = obj$nVessels,     
                          price  = obj$price,           grossSurplus  = obj$grossSurplus,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                          vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                          list(fleet = obj$fleet, year = obj$year, scenario = obj$scenario), 
                          quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
    
      res <- cbind(res[,1:3], 
                 data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                 data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                 data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                 data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]),
                 data.frame(res[,20]), data.frame(res[,21]), data.frame(res[,22]))
                 
      nms1  <- paste('capacity_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('costs_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms8  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms9 <- paste('grossValue_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms10 <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms11 <- paste('fep_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms12 <- paste('netProfit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")     
	  nms13 <- paste('nVessels_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms14 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms15 <- paste('grossSurplus_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms16 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms17 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms18 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms19 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
   
      names(res)[-c(1:3)] <- unlist(mget(paste('nms', 1:19, sep="")))
    }
    else{
        res <- aggregate(list(capacity = obj$capacity,      catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                              discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                              grossValue  = obj$grossValue,         landings  = obj$landings,  fep  = obj$fep, netProfit = obj$netProfit,  nVessels  = obj$nVessels,     
                              price  = obj$price,           grossSurplus  = obj$grossSurplus,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                              vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                         list(fleet = obj$fleet, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
        
        res <- cbind(res[,1:4], 
                     data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),  data.frame(res[,8]),
                     data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),  data.frame(res[,12]),
                     data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]), data.frame(res[,16]),
                     data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]), data.frame(res[,20]),
                     data.frame(res[,21]), data.frame(res[,22]), data.frame(res[,23]))
        
        nms1  <- paste('capacity_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms3  <- paste('costs_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms4  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms5  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms6  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms7  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms8  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms9  <- paste('grossValue_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms10 <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms11 <- paste('fep_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms12 <- paste('netProfit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms13 <- paste('nVessels_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms14 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms15 <- paste('grossSurplus_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms16 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms17 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms18 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms19 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        
        names(res)[-c(1:4)] <- unlist(mget(paste('nms', 1:19, sep="")))
      
    }
  }
  
  return(res)
}


#-------------------------------------------------------------------------------
# revenue_flbeia(fleet, years)
#-------------------------------------------------------------------------------

#' Economic summary functions.
#' 
#' These functions provide summary results of costs, prices and revenues. Provided data can be dessagregated by fleet or by metier depending on the selected function.
#' 
#' @param fleet An element of FLfleets object.
#' @param stock An FLStock object.
#' @param flnm Names of the fleets.
#' @inheritParams FLBEIA
#'   
#' 
#' @details
#'  
#'\itemize{
#'       \item{revenue_flbeia}{ computes the revenue by fleet and metier. The revenue is computed as
#'        landings (weight) multiplied by the price.}
#'       \item{costs_flbeia}{ computes total costs as the sum of fixed and variable costs.}
#'       \item{totvcost_flbeia}{ computes the variable costs including crew share costs .}
#'       \item{totfcost_flbeia}{ computes the total costs by vessel.}
#'       \item{price_flbeia}{ computes the price by stock.} 
#'      }     
#'         


#' @rdname revenue_flbeia
revenue_flbeia <- function(fleet){
    
    sts <- catchNames(fleet)
    mts <- names(fleet@metiers)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    D <- dim(res)
    
    for(mt in mts){
        m <- fleet@metiers[[mt]]
        for(st in sts){
            if(!(st %in% catchNames(m))) next
            dat <- m@catches[[st]]
            res <- res + FLQuant(apply(dat@landings.n*dat@landings.wt*dat@price, c(2,4,6),sum, na.rm=TRUE), dim = D)
        }
    }
    return(res)               
}

#-------------------------------------------------------------------------------
# costs_flbeia(fleet, years)
#-------------------------------------------------------------------------------

#' @rdname revenue_flbeia
#' @aliases costs_flbeia
#' @param covars List of FLQuants with information on covariates.
costs_flbeia <- function(fleet, covars, flnm = NULL){
    
    res <- totvcost_flbeia(fleet) + totfcost_flbeia(fleet, covars, flnm)
    
    return(res)               
}

#-------------------------------------------------------------------------------
# totvcost_flbeia(fleet, years)
#-------------------------------------------------------------------------------
#' @rdname revenue_flbeia
#' @aliases totvcost_flbeia
totvcost_flbeia <- function(fleet){
    
    mts <- names(fleet@metiers)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        res <- res + fleet@metiers[[mt]]@vcost*fleet@effort*fleet@metiers[[mt]]@effshare
    }
    Rev <- revenue_flbeia(fleet)*fleet@crewshare
    
    units(res) <- units(Rev)
    
    res <- res + Rev
    
    return(res)               
}

#-------------------------------------------------------------------------------
# totvcost_flbeia(fleet, years)
#-------------------------------------------------------------------------------
#' @rdname revenue_flbeia
#' @aliases totfcost_flbeia
totfcost_flbeia <- function(fleet, covars, flnm = NULL){
     if(is.null(flnm)) flnm <- 1
     return(fleet@fcost*covars[["NumbVessels"]][flnm, ])            
}


#------------------------------------------------------------------------------#
# fltStkSum :: data.frame[scenario, year, season, stock, fleet, iter, ||,|| 
#        landings, discards, catch, discRat, price, tacshare, quota, quotaUpt] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases fltStkSum
fltStkSum <- function(obj, flnms = names(obj$fleets), stknms = catchNames(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
    
  fleets <- obj$fleets
  advice <- obj$advice
  
#  fleets <- lapply(fleets, setUnitsNA)
 
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')
  
    if(flnms[1] == 'all') flnms <- names(fleets)
    if(stknms[1] == 'all') stknms <- catchNames(fleets)
     
    Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years,])
    
    res <- NULL
    
    
   if(byyear == FALSE){                                
    for(f in flnms){
        
        fl   <- fleets[[f]]

        stfl <- catchNames(fl)        
        sts   <- stknms[stknms %in% stfl]
        
        n <- prod(Dim)*length(sts)
        
        dff <- data.frame(year = rep(years, prod(Dim[2:3])*length(sts)), 
                    season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(sts)), 
                    fleet = rep(f, n), 
                    stock = rep(sts, each = prod(Dim)),
                    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(sts)),  
                    landings = numeric(n), 
                    discards = numeric(n),
                    catch    = numeric(n),
                    discRat  = numeric(n),
                    price    = numeric(n),
                    tacshare = numeric(n),
                    quota    = numeric(n),
                    quotaUpt = numeric(n),
                    stringsAsFactors = FALSE)
        
        k <- 1
        
        for(st in sts){
            
            dff[k:(prod(Dim) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,4,6), sum)[,years])    
            dff[k:(prod(Dim) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,4,6), sum)[,years]) 
            dff[k:(prod(Dim) + k-1),'catch']    <- dff[k:(prod(Dim) + k-1),'discards'] + dff[k:(prod(Dim) + k-1),'landings']
            dff[k:(prod(Dim) + k-1),'discRat']  <- dff[k:(prod(Dim) + k-1),'discards']/dff[k:(prod(Dim) + k-1),'catch']
            dff[k:(prod(Dim) + k-1),'price']    <- c(price_flbeia(fl, st)[,years])
            dff[k:(prod(Dim) + k-1),'tacshare'] <- c((advice$quota.share[[st]][f,])[,years])
            dff[k:(prod(Dim) + k-1),'quota']    <- c((advice$TAC[st,]*advice$quota.share[[st]][f,])[,years])
            dff[k:(prod(Dim) + k-1),'quotaUpt'] <- dff[k:(prod(Dim) + k-1),'catch']/dff[k:(prod(Dim) + k-1),'quota']
            
            k <- k + prod(Dim)     
        }
        res <- rbind(res, dff)
    }}
    else{
      for(f in flnms){
        
        fl   <- fleets[[f]]
        
        stfl <- catchNames(fl)        
        sts  <- stknms[stknms %in% stfl]
        
        n <- prod(Dim[-2])*length(sts)
        
        dff <- data.frame(year = rep(years, prod(Dim[3])*length(sts)), 
                          fleet = rep(f, n), 
                          stock = rep(sts, each = prod(Dim[-2])),
                          iter = rep(rep(1:Dim[3], each = prod(Dim[1])), length(sts)),  
                          landings = numeric(n), 
                          discards = numeric(n),
                          catch    = numeric(n),
                          discRat  = numeric(n),
                          price    = numeric(n),
                          tacshare = numeric(n),
                          quota    = numeric(n),
                          quotaUpt = numeric(n),
                          stringsAsFactors = FALSE)
        
        k <- 1
        
        for(st in sts){
          
          dff[k:(prod(Dim[-2]) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,6), sum)[,years])    
          dff[k:(prod(Dim[-2]) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,6), sum)[,years]) 
          dff[k:(prod(Dim[-2]) + k-1),'catch']    <- dff[k:(prod(Dim[-2]) + k-1),'discards'] + dff[k:(prod(Dim[-2]) + k-1),'landings']
          dff[k:(prod(Dim[-2]) + k-1),'discRat']  <- dff[k:(prod(Dim[-2]) + k-1),'discards']/dff[k:(prod(Dim[-2]) + k-1),'catch']
          dff[k:(prod(Dim[-2]) + k-1),'price']    <- c(seasonMeans(price_flbeia(fl, st)[,years]*quantSums(unitSums(landWStock.f(fl, st)[,years])))/seasonSums(unitSums(quantSums(landWStock.f(fl, st)[,years]))))
          dff[k:(prod(Dim[-2]) + k-1),'tacshare']    <- c((advice$quota.share[[st]][f,])[,years])
          dff[k:(prod(Dim[-2]) + k-1),'quota']    <- c((advice$TAC[st,]*advice$quota.share[[st]][f,])[,years])
          dff[k:(prod(Dim[-2]) + k-1),'quotaUpt'] <- dff[k:(prod(Dim[-2]) + k-1),'catch']/dff[k:(prod(Dim[-2]) + k-1),'quota']
          
          k <- k + prod(Dim[-2])     
        }
        res <- rbind(res, dff)
      }
    }
    
    if(long == TRUE){ # transform res into long format
      r1 <- ifelse(byyear == TRUE, 5,6)
      r2 <- ifelse(byyear == TRUE, 12,13)
      
      names(res)[r1:r2] <- paste('indicator',names(res)[r1:r2], sep = "_")
      res <- reshape(res, direction = 'long', varying = r1:r2, sep = "_")[,1:(r1+1)]
      rownames(res) <- 1:dim(res)[1]
      names(res)[r1:(r1+1)] <- c('indicator', 'value') 
    }
    
    res <- cbind(scenario = scenario, res)
  
    return(res)
}                               

# fltStkSumQ 
#' @rdname bioSum
#' @aliases fltStkSumQ
fltStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario, obj, quantile, prob = prob,na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[6:(6+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario + season, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:6], data.frame(res[,7]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[7:(7+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(catch = obj$catch,  discards = obj$discards, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price,  quota = obj$quota,       quotaUpt = obj$quotaUpt), 
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, scenario = obj$scenario), 
                       quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:4], 
                   data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),  data.frame(res[,8]),
                   data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('quota_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
     
      names(res)[-c(1:4)] <- unlist(mget(paste('nms', 1:7, sep="")))
    }
    else{
      res <- aggregate(list(catch = obj$catch, discards = obj$catch, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price, quota = obj$quota,    quotaUpt = obj$quotaUpt),                   
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:5], 
                   data.frame(res[,6]),  data.frame(res[,7]),   data.frame(res[,8]),  data.frame(res[,9]),
                   data.frame(res[,10]),  data.frame(res[,11]),  data.frame(res[,12]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('quota_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:5)] <-  unlist(mget(paste('nms', 1:7, sep="")))
      
    }
  }
  
  return(res)
}



#-------------------------------------------------------------------------------
# price_flbeia(fleet, years) (mean price in a fleet)
#-------------------------------------------------------------------------------
#' @rdname revenue_flbeia
#' @aliases price_flbeia
price_flbeia <- function(fleet, stock){

    mts <- names(fleet@metiers)
    
    totL <- apply(landWStock.f(fleet, stock), c(2,4,6), sum)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        m <- fleet@metiers[[mt]]
        if(!(stock %in% catchNames(m))) next
        dat <- m@catches[[stock]]
        res <- res + apply(dat@landings.n*dat@landings.wt*dat@price, c(2,4,6),sum, na.rm=TRUE)
    }
    
    res <- res/totL
    
    return(res)                
}


#------------------------------------------------------------------------------#
# mtStkSum data.frame[scenario, year, season, fleet, metier, stock, iter ||,|| 
#        catch, discards, discRat, landings, price] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases mtStkSum
mtStkSum <- function(obj, flnms = names(obj$fleets), stknms = catchNames(obj$fleets), 
                     years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
    
    
  fleets <- obj$fleets
  advice <- obj$advice
  
 # fleets <- lapply(fleets, setUnitsNA)
  
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')
  
  if(flnms[1] == 'all') flnms <- names(fleets)
  if(stknms[1] == 'all') stknms <- catchNames(fleets)
  
    if(flnms[1] == 'all') flnms <- names(fleets)
    if(stknms[1] == 'all') stknms <- catchNames(fleets)
     
    Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years,])

    res <- NULL
    
    if(byyear == FALSE){       
    for(f in flnms){
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        for(m in mts){
            mt   <- fl@metiers[[m]]
            stmt <- catchNames(mt)        
            sts  <- stknms[stknms %in% stmt]

            n <- prod(Dim)*length(sts)
        
            dfm <-  data.frame(year = rep(years, prod(Dim[2:3])*length(sts)), 
                        season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(sts)), 
                        fleet = rep(f, n), 
                        metier = rep(m, n),
                        stock = rep(sts, each = prod(Dim)),
                        iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(sts)),  
                        catch = numeric(n), 
                        discards  = numeric(n),
                        discRat  = numeric(n),
                        landings = numeric(n), 
                        price = numeric(n),
                        stringsAsFactors = FALSE)
            k <- 1
            
            for(ss in sts){
                cc <- mt@catches[[ss]]
                dfm[k:(k+prod(Dim)-1),'landings'] <- c(apply(cc@landings[,years,], c(2,4,6), sum,  na.rm=TRUE))
                dfm[k:(k+prod(Dim)-1),'discards'] <- c(apply(cc@discards[,years,], c(2,4,6), sum,  na.rm=TRUE))
                dfm[k:(k+prod(Dim)-1),'catch'] <- dfm[k:(k+prod(Dim)-1),'discards'] + dfm[k:(k+prod(Dim)-1),'landings']
                dfm[k:(k+prod(Dim)-1),'discRat'] <-  dfm[k:(k+prod(Dim)-1),'discards']/dfm[k:(k+prod(Dim)-1),'catch']
                revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,4,6), sum, na.rm=TRUE)[,years,]
                dfm[k:(k+prod(Dim)-1),'price']  <- c(revst)/dfm[k:(k+prod(Dim)-1),'landings']  
                k <- k + prod(Dim)
            }
            res <- rbind(res, dfm) 
        }  
    }}
    else {
      for(f in flnms){
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        for(m in mts){
          mt   <- fl@metiers[[m]]
          stmt <- catchNames(mt)        
          sts  <- stknms[stknms %in% stmt]
          
          n <- prod(Dim[-2])*length(sts)
          
          dfm <-  data.frame(year = rep(years, prod(Dim[3])*length(sts)), 
                            fleet = rep(f, n), 
                             metier = rep(m, n),
                             stock = rep(sts, each = prod(Dim[-2])),
                             iter = rep(rep(1:Dim[3], each = prod(Dim[1])), length(sts)),  
                             catch = numeric(n), 
                             discards  = numeric(n),
                             discRat  = numeric(n),
                             landings = numeric(n), 
                             price = numeric(n),
                            stringsAsFactors = FALSE)
          k <- 1
          
          for(ss in sts){
            cc <- mt@catches[[ss]]
            dfm[k:(k+prod(Dim[-2])-1),'landings'] <- c(apply(cc@landings[,years,], c(2,6), sum, na.rm=TRUE))
            dfm[k:(k+prod(Dim[-2])-1),'discards'] <- c(apply(cc@discards[,years,], c(2,6), sum, na.rm=TRUE))
            dfm[k:(k+prod(Dim[-2])-1),'catch'] <- dfm[k:(k+prod(Dim[-2])-1),'discards'] + dfm[k:(k+prod(Dim[-2])-1),'landings']
            dfm[k:(k+prod(Dim[-2])-1),'discRat'] <-  dfm[k:(k+prod(Dim[-2])-1),'discards']/dfm[k:(k+prod(Dim[-2])-1),'catch']
            revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,6), sum, na.rm=TRUE)[,years,]
            dfm[k:(k+prod(Dim[-2])-1),'price']  <- c(revst)/dfm[k:(k+prod(Dim[-2])-1),'landings']  
            k <- k + prod(Dim[-2])
          }
          res <- rbind(res, dfm) 
        }  
      }
      
    }
 
    if(long == TRUE){ # transform res into long format
      r1 <- ifelse(byyear == TRUE, 6,7)
      r2 <- ifelse(byyear == TRUE, 10,11)
      
      names(res)[r1:r2] <- paste('indicator',names(res)[r1:r2], sep = "_")
      res <- reshape(res, direction = 'long', varying = r1:r2, sep = "_")[,1:(r1+1)]
      rownames(res) <- 1:dim(res)[1]
      names(res)[r1:(r1+1)] <- c('indicator', 'value') 
    }
    
    res <- cbind(scenario = scenario, res)
    
    return(res)

}


# mtStkSumQ 
#' @rdname bioSum
#' @aliases mtStkSumQ
mtStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + metier + stock + indicator + year + scenario, obj, quantile, prob = prob,na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:6], data.frame(res[,7]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[7:(7+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + metier + stock + indicator + year + scenario + season, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:7], data.frame(res[,8]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[8:(8+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(catch = obj$catch,  discards = obj$discards, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price), 
                       list(fleet = obj$fleet, metier = obj$metier, stock = obj$stock, year = obj$year, scenario = obj$scenario), 
                       quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:5], 
                   data.frame(res[,6]),  data.frame(res[,7]),  data.frame(res[,8]),  
                   data.frame(res[,9]),  data.frame(res[,10]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
      names(res)[-c(1:5)] <- unlist(mget(paste('nms', 1:5, sep="")))
    }
    else{
      res <- aggregate(list(catch = obj$catch, discards = obj$discards, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price),                   
                       list(fleet = obj$fleet, metier = obj$metier, stock = obj$stock, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:6], 
                   data.frame(res[,7]),  data.frame(res[,8]),   data.frame(res[,9]),  data.frame(res[,10]),
                   data.frame(res[,11]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
     
      names(res)[-c(1:6)] <-  unlist(mget(paste('nms', 1:5, sep="")))
      
    }
  }
  
  return(res)
}


#------------------------------------------------------------------------------#
# mtSum data.frame[scenario, year, season, fleet, metier, iter ||,|| 
#        effshare, effort, grossValue, vcost] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases mtSum
mtSum <- function(obj, flnms = names(obj$fleets),
                     years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
  
  if(flnms[1] == 'all') flnms <- names(obj$fleets)
  
  fleets <- obj$fleets
  
 # fleets <- lapply(fleets, setUnitsNA)
  
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')
  
  Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
  Dimnm <- dimnames(fleets[[1]]@effort[,years,])
  
  res <- NULL
  
  if(byyear == FALSE){  
  for(f in flnms){
    fl <- fleets[[f]]
    mts <- names(fl@metiers)
    n <- prod(Dim)*length(mts)
    
    dff <-  data.frame(year = rep(years, prod(Dim[2:3])*length(mts)), 
                       season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(mts)), 
                       fleet = rep(f, n), 
                       metier = rep(mts, each = prod(Dim)),
                       iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(mts)),  
                       effshare = numeric(n), 
                       effort = numeric(n),
                       grossValue = numeric(n), 
                       vcost = numeric(n),
                       stringsAsFactors = FALSE)
    k <- 1
    for(m in mts){
      mt <- fl@metiers[[m]]
      dff[k:(k+prod(Dim)-1),'effort']   <- c((fl@effort*mt@effshare)[,years,])
      dff[k:(k+prod(Dim)-1),'effshare'] <- c(mt@effshare[,years,])
      dff[k:(k+prod(Dim)-1),'vcost']    <- c((fl@effort*mt@effshare*mt@vcost)[,years,])
       dff[k:(k+prod(Dim)-1),'grossValue'] <- c(Reduce('+', lapply(mt@catches, function(x) unitSums(quantSums(x@landings.n*x@price))[,years])))

      k <- k + prod(Dim)
    }
    
    res <- rbind(res, dff)      
  }}
  else{
    for(f in flnms){
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      n <- prod(Dim[-2])*length(mts)
      
      dff <-  data.frame(year = rep(years, Dim[3]*length(mts)), 
                         fleet = rep(f, n), 
                         metier = rep(mts, each = prod(Dim[-2])),
                         iter = rep(rep(1:Dim[3], each = prod(Dim[1])), length(mts)),  
                         effshare = numeric(n), 
                         effort = numeric(n),
                         grossValue = numeric(n), 
                         vcost = numeric(n),
                         stringsAsFactors = FALSE)
      k <- 1
      for(m in mts){
        mt <- fl@metiers[[m]]
        dff[k:(k+prod(Dim[-2])-1),'effort']   <- c(seasonSums((fl@effort*mt@effshare)[,years,]))
        dff[k:(k+prod(Dim[-2])-1),'effshare'] <- c(seasonSums(mt@effshare[,years,]))
        dff[k:(k+prod(Dim[-2])-1),'vcost']    <- c(seasonSums(fl@effort*mt@effshare*mt@vcost)[,years,])
        dff[k:(k+prod(Dim[-2])-1),'grossValue'] <- c(Reduce('+', lapply(mt@catches, function(x) seasonSums(unitSums(quantSums(x@landings.n*x@price)))[,years])))
        
        k <- k + prod(Dim[-2])
      }
      
      res <- rbind(res, dff)   
      
    }
  }
    
  if(long == TRUE){
     # transform res into long format
      r1 <- ifelse(byyear == TRUE, 5,6)
      r2 <- ifelse(byyear == TRUE, 8,9)
      
      names(res)[r1:r2] <- paste('indicator',names(res)[r1:r2], sep = "_")
      res <- reshape(res, direction = 'long', varying = r1:r2, sep = "_")[,1:(r1+1)]
      rownames(res) <- 1:dim(res)[1]
      names(res)[r1:(r1+1)] <- c('indicator', 'value') 
    }
    
    res <- cbind(scenario = scenario, res)
    

  return(res)
}


# mtStkSumQ 
#' @rdname bioSum
#' @aliases mtSumQ
mtSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 9){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + metier + indicator + year + scenario, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[6:(6+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + metier + indicator + year + scenario + season, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:6], data.frame(res[,7]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[7:(7+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(effort = obj$effort,  effshare = obj$effshare, vcost = obj$vcost, grossValue = obj$grossValue,       
                            vcost = obj$vcost), 
                       list(fleet = obj$fleet, metier = obj$metier, year = obj$year, scenario = obj$scenario), 
                       quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:4], 
                   data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),  
                   data.frame(res[,6]))
      
      nms1  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('effshare_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('grossValue_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
       
      names(res)[-c(1:4)] <- unlist(mget(paste('nms', 1:4, sep="")))
    }
    else{
      res <- aggregate(list(effort = obj$effort,  effshare = obj$effshare, vcost = obj$vcost, grossValue = obj$grossValue,       
                            vcost = obj$vcost),                    
                       list(fleet = obj$fleet, metier = obj$metier,  year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:5], 
                   data.frame(res[,6]),  data.frame(res[,7]),   data.frame(res[,8]),  data.frame(res[,9]))
                   
      
      nms1  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('effshare_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('grossValue_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:5)] <-  unlist(mget(paste('nms', 1:4, sep="")))
      
    }
  }
  
  return(res)
}

#------------------------------------------------------------------------------#
# advSum :: data.frame[scenario, year, stock, iter ||,|| 
#        catch, discards, discRat, landings, quotaUpt, tac] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases advSum

advSum <- function(obj, stknms = 'all', years = dimnames(obj$biols[[1]]@n)$year, long = FALSE, scenario = 'bc'){
  
  if(stknms == 'all') stknms <- names(obj$biols)  
  
  x1 <- Reduce(rbind, lapply(stknms, function(x)  cbind(stock = x, 
                                                               array2df(apply(catchWStock(obj$fleets, x), c(2,6), sum), label.x = 'catch')[,c('year', 'iter', 'catch')])))
  
  x2 <- Reduce(rbind, lapply(stknms, function(x)  cbind(stock = x, 
                                                               array2df(apply(landWStock(obj$fleets, x), c(2,6), sum), label.x = 'landings')[,c('year', 'iter', 'landings')])))
  
  x3 <- Reduce(rbind, lapply(stknms, function(x)  cbind(stock = x, 
                                                               array2df(apply(discWStock(obj$fleets, x), c(2,6), sum), label.x = 'discards')[,c('year', 'iter', 'discards')])))
  
  res <- as.tbl(cbind(x1,discards = x3[,4], landings = x3[,4]))
                
  x4 <- as.tbl(array2df(obj$advice$TAC, label.x = 'tac')[,c('stock', 'year', 'iter', 'tac')])
                
  res <- full_join(res, x4, by = c('stock', 'year', 'iter'))
                
  # Wide format 
  res <- res %>%  dplyr::group_by(stock, year, iter) %>% mutate(quotaUpt = catch/tac, discRat = discards/catch, scenario = scenario) 
                
  # reshaping this to the long format
   if(long == TRUE) res <- res %>% gather(key=indicator, value=value, catch, discards, discRat,landings, quotaUpt, tac)
                
 return(res)
}

#' @rdname advSum
#' @aliases advSumQ
advSumQ <- function(obj,  probs = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] <= 7){ # the object is in long format
    
    res <- obj %>% dplyr::group_by(stock, year, indicator)  %>% 
      dplyr::summarise(qlow = quantile(value, probs=probs[1], na.rm = TRUE), 
                qmed = quantile(value, probs=probs[2], na.rm = TRUE), 
                qupp = quantile(value, probs=probs[3], na.rm = TRUE))
    
  }
  else{
    
    res1 <- obj %>% dplyr::group_by(stock, year) %>%  
      summarise_at(c("catch",    "discards", "discRat",  "landings", "quotaUpt", "tac"),
                   .funs =  list(qlow = quantile),probs= probs[1], na.rm=T)
    
    res2 <- obj %>% dplyr::group_by(stock, year) %>%  
      summarise_at(c("catch",    "discards", "discRat",  "landings", "quotaUpt", "tac"),
                   .funs =  list(qmed = quantile),probs=probs[2], na.rm=T)
    
    res3 <- obj %>% dplyr::group_by(stock, year) %>%  
      summarise_at(c("catch",    "discards", "discRat",  "landings", "quotaUpt", "tac"),
                   .funs =  list(qupp = quantile),probs=probs[3], na.rm=T)
    
    res <- bind_cols(res1, res2[,-(1:2)], res3[,-(1:2)])
    
  }
  
  return(res)
}

#----------------------------------------------------------------------
# riskSum(obj, stocks, fleets, years, long)
# Bpa = a named vector with the precautionary biomass per stock.
# Blim = a named vector with the limit biomass per stock.
# Prflim = a named vector with the limit profit per fleet.
#----------------------------------------------------------------------
#' @rdname bioSum
#' @aliases riskSum
riskSum <- function(obj, stknms = names(obj$biols), Bpa, Blim, Prflim, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], scenario = 'bc'){
  
  # biols
  
  if (stknms == 'all') { 
    stknms <- names(obj$biols)
  } else if (sum(!stknms %in% names(obj$biols))>0) {
    stop(paste("'stknms' values should be in the following list: ", paste(names(obj$biols), collapse = ", "), sep=''))
  }
  if (sum(!names(Bpa) %in% stknms) + sum(!names(Blim) %in% stknms)>0) {
    stop(paste("Check names for 'Bpa' and 'Blim'. Values should be in the following list: ", paste(stknms, collapse = ", "), sep=''))
  }
  
  bioS <- bioSum(obj, stknms = stknms, years = years, long = FALSE, scenario = scenario)
  
  bioS <- bioS %>% dplyr::group_by(scenario, year, stock, iter) %>% 
    mutate(Bpa = Bpa[stock], Blim = Blim[stock], risk.pa = as.numeric(ssb<Bpa), risk.lim = as.numeric(ssb<Blim))
  
  bioS.pa <- bioS %>% dplyr::group_by(year, stock, scenario) %>% 
    dplyr::summarise(indicator="pBpa", value = sum(risk.pa)/length(risk.pa))
  
  bioS.lim <- bioS %>% dplyr::group_by(year, stock, scenario) %>% 
    dplyr::summarise(indicator = "pBlim", value = sum(risk.lim)/length(risk.lim))
  
  outbio <- bind_rows( bioS.lim, bioS.pa) %>% dplyr::rename(unit=stock)
  
  # fleets
  
  if (flnms == 'all') { 
    flnms <- names(obj$fleets)
  } else if (sum(!flnms %in% names(obj$fleets))>0) {
    stop(paste("'flnms' values should be in the following list: ", paste(names(obj$fleets), collapse = ", "), sep=''))
  }
  if (sum(!names(Prflim) %in% flnms)>0) {
    stop(paste("Check names for 'Prflim'. Values should be in the following list: ", paste(flnms, collapse = ", "), sep=''))
  }
  
  flS <- fltSum(obj, years = years, flnms = flnms, long = FALSE, scenario = scenario)
  
  flS <- flS %>% dplyr::group_by(scenario, year, fleet, iter) %>% 
    mutate(refp = Prflim[fleet], risk = as.numeric(grossSurplus < refp))
  
  outfl <- flS %>% dplyr::group_by(year, fleet, scenario) %>% 
    dplyr::summarise(indicator = "pPrflim", value = sum(risk)/length(risk)) %>% dplyr::rename(unit=fleet)
  
  # all combined
  
  res <- bind_rows( outbio, outfl)
  
  return(res)
}


#----------------------------------------------------------------------
# npv(obj, years, flnms)
#----------------------------------------------------------------------
#' @rdname bioSum
#' @aliases npv
npv <- function(obj, discF = 0.05, y0, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], scenario = 'bc'){
  
  flS <- fltSum(obj, years = years, flnms = flnms, long = FALSE)
  
  y0 <- as.numeric(y0)
  
  flS <- cbind(flS, discount= (1+discF)^(as.numeric(as.character(flS$year))-y0))
  flS <- cbind(flS, discProf = flS$grossSurplus/(flS$discount))
  
  flS <- subset(flS, year %in% c(years))
  
  res <- aggregate(discProf ~ fleet + iter, data = flS, FUN = sum)
  
  names(res)[3] <- 'npv'
  
  res <- cbind(scenario = scenario, res)
  
  return(res)
  
}  
    

#' @rdname bioSum
#' @aliases npvQ
npvQ <- function(obj, prob = c(0.05,0.5,0.95)){
  
  res <- aggregate(npv ~  fleet + scenario, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
  res <- cbind(res[,1:2], data.frame(res[,3]))
  
  nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
  
  names(res)[3:5] <- nms
  
  
  return(res)
  
}   


#----------------------------------------------------------------------
# vesselSum
#----------------------------------------------------------------------
#' @rdname bioSum
#' @aliases vesselSum
vesselSum <- function(obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = TRUE, scenario = 'bc'){
  
  flS <- fltSum(obj, flnms = flnms, years = years, byyear = byyear, long = long, scenario = scenario)
  
  ids <- c("catch", "costs", "discards","discRat","effort",  "fcosts","gva","grossValue", "landings",     
           "fep",  "netProfit", "price", "grossSurplus", "quotaUpt", "salaries", "vcosts")
  
  if(byyear == TRUE){
  if(long == FALSE){
    for(col in ids){ #5:22
      flS[,col] <- flS[,col]/flS[,'nVessels']
    }
    res <- flS[,c(names(flS)[1:4], ids, 'profitability')] # res <- flS[,c(1:4,6:15,17:23)]
  }
  else{
   # ids <- c("catch", "costs", "discards","discRat","effort",  "fcosts","gva","grossValue", "landings",     
   #    "fep",  "netProfit", "price", "grossSurplus", "quotaUpt", "salaries", "vcosts")
   for(id in ids){
     flS[flS$indicator == id, 'value'] <- flS[flS$indicator == id, 'value']/flS[flS$indicator == 'nVessels', 'value']
   }
   res <- subset(flS, indicator %in% c(ids, 'profitability'))
  }}
   else{
     if(long == FALSE){
       for(col in ids){ #6:23
         flS[,col] <- flS[,col]/flS[,'nVessels']
       }
       res <- flS[,c(names(flS)[1:5], ids, 'profitability')] # res <- flS[,c(1:5,7:16,18:24)]
     }
     else{
       # ids <- c("catch", "costs", "discards","discRat","effort",  "fcosts","gva","grossValue", "landings",     
       #          "fep", "netProfit",  "price", "grossSurplus", "quotaUpt", "salaries", "vcosts")
       for(id in ids){
         flS[flS$indicator == id, 'value'] <- flS[flS$indicator == id, 'value']/flS[flS$indicator == 'nVessels', 'value']
       }
       res <- subset(flS, indicator %in% c(ids, 'profitability'))
   }
   
  }
     return(res) 
               
}

#' @rdname bioSum
#' @aliases vesselSumQ
vesselSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){

if(dim(obj)[2] < 10){ # the object is in long format
  
  if(!('season' %in% names(obj))){
    res <- aggregate(value ~ fleet + indicator + year + scenario, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
    names(res)[5:(5+length(prob)-1)] <- nms
  }
  else{
    res <- aggregate(value ~ fleet + indicator + year + scenario + season, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
    res <- cbind(res[,1:5], data.frame(res[,6]))
    
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[6:(6+length(prob)-1)] <- nms
  }
}
else{
  
  if(!('season' %in% names(obj))){
    res <- aggregate(list(    catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                          discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                          grossValue  = obj$grossValue,         landings  = obj$landings,  fep  = obj$fep,  netProfit = obj$netProfit,  
                          price  = obj$price,           grossSurplus  = obj$grossSurplus,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                          vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                     list(fleet = obj$fleet, year = obj$year, scenario = obj$scenario), 
                     quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
    
    res <- cbind(res[,1:3], 
                 data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                 data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                 data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                 data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]),  
                 data.frame(res[,20]))
    
    nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms2  <- paste('costs_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms3  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms4  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms5  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms6  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms7  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms8  <- paste('grossValue_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms9  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms10 <- paste('fep_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms11 <- paste('netProfit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms12 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms13 <- paste('grossSurplus_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms14 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms15 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms16 <- paste('vcosts_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms17 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
    names(res)[-c(1:3)] <- unlist(mget(paste('nms', c(1:17), sep="")))
  }
  else{
    res <- aggregate(list(    catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                          discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                          grossValue  = obj$grossValue,         landings  = obj$landings,  fep  = obj$fep, netProfit = obj$netProfit,   
                          price  = obj$price,           grossSurplus  = obj$grossSurplus,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                          vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                     list(fleet = obj$fleet, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
    
    res <- cbind(res[,1:4], 
                 data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]), data.frame(res[,8]),  
                 data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]), data.frame(res[,12]), 
                 data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]), data.frame(res[,16]), 
                 data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]), data.frame(res[,20]), 
                 data.frame(res[,21]))
    
    nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms2  <- paste('costs_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms3  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms4  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms5  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms6  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms7  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms8  <- paste('grossValue_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms9  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms10 <- paste('fep_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms11 <- paste('netProfit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms12 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms13 <- paste('grossSurplus_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms14 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms15 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms16 <- paste('vcosts_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms17 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
    names(res)[-c(1:4)] <- unlist(mget(paste('nms',  c(1:17), sep="")))
    
  }
}

return(res)
}


#----------------------------------------------------------------------
# vesselStkSum
#----------------------------------------------------------------------
#' @rdname bioSum
#' @aliases vesselStkSum
vesselStkSum <- function(obj, flnms = names(obj$fleets), stknms = catchNames(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
  
  fleets <- obj$fleets
  covars <- obj$covars
  
  if (flnms[1] == "all") flnms <- names(fleets)
  
  Dim <- dim(fleets[[1]]@effort[, years, ])[c(2, 4, 6)]
  Dimnm <- dimnames(fleets[[1]]@effort[, years, ])
  n <- prod(Dim) * length(flnms)
  
  if(is.null(covars$Depreciation)){ 
    covars$Depreciation <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, years, ])[2:6]))
    dimnames(covars$Depreciation) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, years, ])[2:6])
  }
  if(is.null(covars$Salaries)){ 
    covars$Salaries <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, years, ])[2:6]))
    dimnames(covars$Salaries) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, years, ])[2:6])
  }
  if(is.null(covars$MaxDays)){ 
    covars$MaxDays <- FLQuant(365/dim(fleets[[1]]@effort[, years, ])[4], dim = c(length(fleets),dim(fleets[[1]]@effort[, years, ])[2:6]))
    dimnames(covars$MaxDays) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, years, ])[2:6])
  }
  if(is.null(covars$NumbVessels)){ 
    covars$NumbVessels <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, years, ])[2:6]))
    dimnames(covars$NumbVessels) <- c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, years, ])[2:6])
  }
  
  flS <- fltStkSum(obj, flnms = flnms, stknms = stknms, years = years, byyear = byyear, long = long, scenario = scenario)
  
  if(byyear == TRUE){
    if(long == FALSE){
      for(col in c(6:8,10:12)){
        for(fl in flnms){
          for(st in stknms){
            flS[flS$fleet == fl & flS$stock == st,col] <- flS[flS$fleet == fl & flS$stock == st,col]/c(seasonMeans(covars[['NumbVessels']][fl,years]))
          }}}
      res <- flS
    }
    else{
      ids <- c("landings", "discards", "catch" ,   "price",  "tacshare",   "quota"   )
      for(id in ids){
        for(fl in flnms){
          for(st in stknms){
            flS[flS$indicator == id & flS$fleet == fl & flS$stock == st, 'value'] <- flS[flS$indicator == id & flS$fleet == fl & flS$stock == st, 'value']/c(seasonMeans(covars[['NumbVessels']][fl,years]))
      }}}
  
      res <- flS
    }}
  if(byyear == FALSE){
    if(long == FALSE){
      for(col in c(6:8,10:12)){
        for(fl in flnms){
          for(st in stknms){
            for(ss in dimnames(fleets[[fl]]@effort)[[4]]){
              flS[flS$fleet == fl & flS$stock == st & flS$season == ss,col] <- flS[flS$fleet == fl & flS$stock == st & flS$season == ss,col]/c((covars[['NumbVessels']][fl,years,,ss]))
          }}}}
      res <- flS
    }
    else{
      ids <- c("landings", "discards", "catch" ,   "price",  "tacshare",   "quota"   )
      for(id in ids){
        for(fl in flnms){
          for(st in stknms){
            for(ss in dimnames(fleets[[fl]]@effort)[[4]]){
              flS[flS$indicator == id & flS$fleet == fl & flS$stock == st & flS$season == ss, 'value'] <- flS[flS$indicator == id & flS$fleet == fl & flS$stock == st & flS$season == ss, 'value']/c(seasonMeans(covars[['NumbVessels']][fl,years,,ss]))
          }}}}
      
      res <- flS
    }}
  
  return(res) 
  
}

# vesselStkSumQ 
#' @rdname bioSum
#' @aliases vesselStkSumQ
vesselStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[6:(6+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario + season, obj, quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      res <- cbind(res[,1:6], data.frame(res[,7]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[7:(7+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(catch = obj$catch,  discards = obj$discards, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price,  quota = obj$quota,       quotaUpt = obj$quotaUpt), 
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, scenario = obj$scenario), 
                       quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:4], 
                   data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),  data.frame(res[,8]),
                   data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('quota_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:4)] <- unlist(mget(paste('nms', 1:7, sep="")))
    }
    else{
      res <- aggregate(list(catch = obj$catch, discards = obj$catch, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price, quota = obj$quota,    quotaUpt = obj$quotaUpt),                   
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.action = na.pass, na.rm=TRUE)
      
      res <- cbind(res[,1:4], 
                   data.frame(res[,5]),  data.frame(res[,6]),   data.frame(res[,7]),  data.frame(res[,8]),
                   data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('quota_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:4)] <-  unlist(mget(paste('nms', 1:7, sep="")))
      
    }
  }
  
  return(res)
}


