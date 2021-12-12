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
#' @description These functions return the biomass (B), fishing mortality (F),  spawning stock biomass (SSB), recruitment (R), catches (C), landings (L) and discards (D) indicators. Also indicators comparing the reference points with the actual values, Bpa, Blim, Btarget, Fpa, Flim and Ftarget, so the biomass indicators are TRUE if the biomass is above them, and the fishing mortality indicators are TRUE if the fishing mortality is below them. ssb2Btarget and f2Ftarget return the ratio between SSB and F and the target reference point.   
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
            m.  <- array(unitMeans(obj$biols[[stk]]@m)[fbar_age,years,,,drop=T], dim = c(length(fbar_age),ny,ns,it), dimnames = Dnms)
            c.  <- array(apply(catchStock(obj$fleets, stk),c(1:2,4,6), sum)[fbar_age,years,,,drop = TRUE], dim = c(length(fbar_age),ny,ns,it), dimnames = Dnms)
        
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
                          aux[a,y,ss,i] <- ifelse(class(xx) == 'try-error', NA, xx)
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
    
    res <- array(0,dim = c(length(stknms), ny,ns,it), dimnames = list(stock = stknms, year = yrnms, season = ssnms, iter = 1:it))
    
    for(stk in stknms){
      if(dim(obj$biols[[stk]]@n)[[1]] > 1) for(ss in 1:ns) {
        uu <- ifelse( dim(obj$biols[[stk]]@n)[[3]] > 1, ss, 1)
        res[stk,,ss,] <- obj$biols[[stk]]@n[1,yrnms,uu,ss,drop=T]
      }
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
#' @param long logical. The data frame should be constructed using long or wide format? Default FALSE.
#' @param byyear logical. The indicators should be provided at season or year level? Default TRUE.
#' @param ssb_season If byyear = TRUE, the season in which ssb will be taken.
#' @param prob a numeric vector with the probabilities used to calculate the quantiles. 
#' @param scenario a character string with the name of the scenario corresponding with obj. Default bc.
#' @param Bpa named numeric vector with one element per stock in stknms. The precautionary approach stock spawning biomass used in riskSum function to calculate biological risk yearly.
#' @param Blim named numeric vector with one element per stock in stknms. The limit stock spawning biomass used in riskSum function to calculate biological risk yearly.
#' @param ProfRef named numeric vector with one element per fleet in flnms. The reference profit level used in riskSum function to calculate economic risk yearly.
#' @param discF Discount rate.
#' @param y0 character. Reference year.
#' @param verbose logical. If TRUE, prints the function steps.
#' @param  brp a data frame with columns stock, iter and one colum per reference point with the value of the biological reference points per stock and iteration. 
#' The used reference points are Bpa, Blim, Btarget, Fpa, Flim and Ftarget. 

#' @examples
#'\dontrun{
#'
#' library(FLBEIA)
#'
#' # Apply the summary functions to the examples runs in FLBEIA help page.
#' # Test the different arguments in summary function.
#' 
#' data(res_flbeia)
#' 
#' 
#' #------------------------------------------------
#' # Example One: One stock, one fleet, one iter.
#' #------------------------------------------------
#' 
#' # Wide format (default)
#' 
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
#' # Long format for a range of years
#' 
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
#' # Wide format with seasonal disaggregation 
#' # (Note: No seasonal disagregation available for adv summaries)
#' 
#' oneRes_bio    <- bioSum(oneRes, byyear = FALSE)
#' oneRes_flt    <- fltSum(oneRes, byyear = FALSE)
#' oneRes_fltStk <- fltStkSum(oneRes, byyear = FALSE)
#' oneRes_mt     <- mtSum(oneRes, byyear = FALSE)
#' oneRes_mtStk  <- mtStkSum(oneRes, byyear = FALSE)
#' oneRes_adv    <- advSum(oneRes) # Advice summary is only by year.
#' 
#' oneRes_bioQ    <- bioSumQ(oneRes_bio)
#' oneRes_fltQ    <- fltSumQ(oneRes_flt)
#' oneRes_fltStkQ <- fltStkSumQ(oneRes_fltStk)
#' oneRes_mtQ     <- mtSumQ(oneRes_mt)
#' oneRes_mtStkQ  <- mtStkSumQ(oneRes_mtStk)
#' oneRes_advQ    <- advSumQ(oneRes_adv)
#' 
#' # Long format and seasonal
#' 
#' oneRes_bio    <- bioSum(oneRes, long = TRUE)
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
#' 
#' # Wide format (default)
#' 
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
#' # Long format for a range of years
#' 
#' oneItRes_bio    <- bioSum(oneItRes, long = TRUE, years = ac(2016:2020))
#' oneItRes_flt    <- fltSum(oneItRes, long = TRUE, years = ac(2016:2020))
#' oneItRes_fltStk <- fltStkSum(oneItRes, long = TRUE, years = ac(2016:2020))
#' oneItRes_mt     <- mtSum(oneItRes, long = TRUE, years = ac(2016:2020))
#' oneItRes_mtStk  <- mtStkSum(oneItRes, long = TRUE, years = ac(2016:2020))
#' oneItRes_adv    <- advSum(oneItRes, long = TRUE, years = ac(2016:2020))
#' 
#'
#' oneItRes_bioQ    <- bioSumQ(oneItRes_bio)
#' oneItRes_fltQ    <- fltSumQ(oneItRes_flt)
#' oneItRes_fltStkQ <- fltStkSumQ(oneItRes_fltStk)
#' oneItRes_mtQ     <- mtSumQ(oneItRes_mt)
#' oneItRes_mtStkQ  <- mtStkSumQ(oneItRes_mtStk)
#' oneItRes_advQ    <- advSumQ(oneItRes_adv)
#' 
#' # Wide format with seasonal disaggregation 
#' # (Note: No seasonal disagregation available for adv summaries)
#' 
#' oneItRes_bio    <- bioSum(oneItRes, byyear = FALSE)
#' oneItRes_flt    <- fltSum(oneItRes, byyear = FALSE)
#' oneItRes_fltStk <- fltStkSum(oneItRes, byyear = FALSE)
#' oneItRes_mt     <- mtSum(oneItRes, byyear = FALSE)
#' oneItRes_mtStk  <- mtStkSum(oneItRes, byyear = FALSE)
#' oneItRes_adv    <- advSum(oneItRes) # Advice summary is only by year.
#' 
#' oneItRes_bioQ    <- bioSumQ(oneItRes_bio)
#' oneItRes_fltQ    <- fltSumQ(oneItRes_flt)
#' oneItRes_fltStkQ <- fltStkSumQ(oneItRes_fltStk)
#' oneItRes_mtQ     <- mtSumQ(oneItRes_mt)
#' oneItRes_mtStkQ  <- mtStkSumQ(oneItRes_mtStk)
#' oneItRes_advQ    <- advSumQ(oneItRes_adv)
#' 
#' # Long format and seasonal
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
#' 
#' #------------------------------------------------
#' # Example Multi: Two stock, two fleet, four iters.
#' #------------------------------------------------
#' 
#' # Wide format (default)
#' 
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
#' # Long format for a range of years
#' 
#' multiRes_bio    <- bioSum(multiRes, long = TRUE, years = ac(2016:2020))
#' multiRes_flt    <- fltSum(multiRes, long = TRUE, years = ac(2016:2020))
#' multiRes_fltStk <- fltStkSum(multiRes, long = TRUE, years = ac(2016:2020))
#' multiRes_mt     <- mtSum(multiRes, long = TRUE, years = ac(2016:2020))
#' multiRes_mtStk  <- mtStkSum(multiRes, long = TRUE, years = ac(2016:2020))
#' multiRes_adv    <- advSum(multiRes, long = TRUE, years = ac(2016:2020))
#' 
#' multiRes_bioQ    <- bioSumQ(multiRes_bio)
#' multiRes_fltQ    <- fltSumQ(multiRes_flt)
#' multiRes_fltStkQ <- fltStkSumQ(multiRes_fltStk)
#' multiRes_mtQ     <- mtSumQ(multiRes_mt)
#' multiRes_mtStkQ  <- mtStkSumQ(multiRes_mtStk)
#' multiRes_advQ    <- advSumQ(multiRes_adv)
#' 
#' # Wide format with seasonal disaggregation 
#' # (Note: No seasonal disagregation available for adv summaries)
#' 
#' multiRes_bio    <- bioSum(multiRes, byyear = FALSE)
#' multiRes_flt    <- fltSum(multiRes, byyear = FALSE)
#' multiRes_fltStk <- fltStkSum(multiRes, byyear = FALSE)
#' multiRes_mt     <- mtSum(multiRes, byyear = FALSE)
#' multiRes_mtStk  <- mtStkSum(multiRes, byyear = FALSE)
#' multiRes_adv    <- advSum(multiRes) # Advice summary is only by year.
#' 
#' multiRes_bioQ    <- bioSumQ(multiRes_bio)
#' multiRes_fltQ    <- fltSumQ(multiRes_flt)
#' multiRes_fltStkQ <- fltStkSumQ(multiRes_fltStk)
#' multiRes_mtQ     <- mtSumQ(multiRes_mt)
#' multiRes_mtStkQ  <- mtStkSumQ(multiRes_mtStk)
#' multiRes_advQ    <- advSumQ(multiRes_adv)
#' 
#' # Long format and seasonal
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
bioSum <- function(obj, stknms = 'all', years = dimnames(obj$biols[[1]]@n)$year, long = FALSE, scenario = 'bc', byyear = TRUE, ssb_season = NULL, brp = NULL){
  
  # For avoiding warnings in R CMD CHECK
  year <- season <- NULL
  
  if(stknms == 'all') stknms <- names(obj$biols)  
  
  if(byyear == TRUE & is.null(ssb_season)) ssb_season <- dimnames(obj$biols[[1]]@n)[[4]][1]
  
  xx <- summary_flbeia(obj) # array: stk x year x season x iter x indicator
  dat <- array2df(xx, label.x="value")
  
  # Wide format 
  res <- dat %>% mutate(scenario = scenario) %>% tidyr::spread(key='indicators', value='value') %>% 
    filter(year %in% years & stock %in% stknms) %>% arrange(year) %>%
    mutate_at(vars(c(2,4)), as.character) %>% mutate_at(vars(c(2,4)), as.numeric) %>%
    dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$season, .data$iter) %>% 
    mutate(catch.iyv = catch/dplyr::lag(catch), land.iyv = landings/dplyr::lag(landings), disc.iyv = discards/dplyr::lag(discards)) 
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "stock"     "year"      "season"    "iter"      "scenario"  "rec":"disc.iyv"
  # character   numeric     character   numeric     character   numeric
  
  res <- res %>% mutate(stock = as.character(stock),
                        season = as.character(season))
  
  # year or seasonal?
  if(byyear == TRUE){
    if(length(unique(res$season)) >1){
      # indicators that are summ up over the seasons
      res1 <- res %>% dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$iter)  %>% summarise_at(c('rec','f','catch', 'landings', 'discards'),'sum')
      # ssb user selects the season
      res2 <- res %>%  filter(season == ssb_season) %>% ungroup() %>% select('stock', 'year', 'iter', 'scenario', 'ssb')
      # biomass the first season
      res3 <- res %>%  filter(season == dimnames(obj$biols[[1]]@n)[[4]][1]) %>% ungroup() %>% select('stock', 'year', 'iter', 'scenario', 'biomass')
      
      res <- full_join(res1, res2, by = c('stock', 'year', 'iter', 'scenario'))
      res <- full_join(res, res3, by = c('stock', 'year', 'iter', 'scenario'))
      res <- res %>% arrange(year) %>%  dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$iter) %>%  mutate(catch.iyv = catch/dplyr::lag(catch), land.iyv = landings/dplyr::lag(landings),
                                                                                                                         disc.iyv = discards/dplyr::lag(discards)) 
    }
    else{
      res <- res %>% ungroup() %>% select(-season)
    }
  }
  
  # If brp is not provided create it with NAs.
  if(is.null(brp)) brp <- as_tibble(cbind(expand.grid(stock = unique(res$stock), iter = unique(res$iter)),
                                          Ftarget = NA, Btarget = NA, Flim = NA, Fpa = NA, Blim = NA, Bpa = NA))
  
  brp <- as_tibble(brp) %>% ungroup() %>% group_by(stock, iter)

  res <- res %>%  ungroup() %>% group_by(stock, iter) %>% left_join(brp)
  res <- res %>% mutate(ssb2Btarget = ssb/Btarget, f2Ftarget = f/Ftarget,
                        lowerBpa  = ifelse(ssb<Bpa, TRUE, FALSE), higherFpa = ifelse(f>Fpa, TRUE, FALSE),
                        lowerBlim = ifelse(ssb<Blim, TRUE, FALSE), higherFlim = ifelse(f>Flim, TRUE, FALSE),
                        lowerBtarget = ifelse(ssb<Btarget, TRUE, FALSE), higherFtarget = ifelse(f>Ftarget, TRUE, FALSE))

  
  # reshaping this to the long format
  if(long == TRUE){
    res <- res %>% gather(key='indicator', value='value', c("biomass", "f", "rec", "ssb", "catch", "landings", "discards", "catch.iyv", "land.iyv", "disc.iyv",
                                                            "higherFtarget", "lowerBtarget", "higherFlim", "higherFpa", "lowerBlim","lowerBpa",  "ssb2Btarget", "f2Ftarget"))
  }
  
  return(res)
}


#' @rdname bioSum
#' @aliases bioSumQ

bioSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  p_names <- paste("q",ifelse(nchar(substr(prob,3, nchar(prob)))==1, 
                              paste(substr(prob,3, nchar(prob)), 0, sep = ""), 
                              substr(prob,3, nchar(prob))), sep = "")
  
  if("indicator" %in% names(obj)){ # the object is in long format
    
    objRP <- obj %>% filter(indicator %in% c('higherFpa', 'higherFlim', 'lowerBpa', 'lowerBlim', 'higherFtarget', 'lowerBtarget'))
    obj   <- obj %>% filter(!(indicator %in% c('higherFpa', 'higherFlim', 'lowerBpa', 'lowerBlim', 'higherFtarget', 'lowerBtarget')))
    
    if('season' %in% names(obj)){
      res <- obj %>% dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$season, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
      
      # BRP indicators
      resRP <- objRP %>% dplyr::group_by(.data$scenario,.data$year,.data$season, .data$stock, .data$indicator) %>%
        dplyr::summarise(value = sum(value, na.rm=T)/dplyr::n()) 
      
      resRP <- bind_cols(resRP[,1:5],NA, resRP[,6], NA)
      names(resRP)[6:8] <- names(res)[6:8]
      
    }
    
    else{
      res <- obj %>% dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$indicator) %>% 
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
      
      # BRP indicators
      resRP <- objRP %>% dplyr::group_by(.data$scenario,.data$year, .data$stock, .data$indicator) %>%
        dplyr::summarise(value = sum(value)/dplyr::n()) 
      
      resRP <- bind_cols(resRP[,1:4],NA, resRP[,5], NA)
      names(resRP)[5:7] <- names(res)[5:7]
    }
    res   <- bind_rows(res, resRP)
    
    res$indicator <- dplyr::recode(res$indicator, higherFpa = 'pFpa', higherFlim = 'pFlim', lowerBpa = 'pBpa', 
                                           lowerBlim = 'pBlim', higherFtarget = 'pFtarget', lowerBtarget = 'pBtarget')
  }
  else{
    
    p_funs <- purrr::map(prob, ~purrr::partial(quantile, prob = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    if('season' %in% names(obj)){
      
      res <- obj %>% dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$season) %>%  
        summarise_at(c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards', 
                       'catch.iyv', 'land.iyv', 'disc.iyv', 'ssb2Btarget', 'f2Ftarget'),
                     .funs =  p_funs)
      # RP indicator
      res <- res %>% mutate(pBlim_q95 = NA, pBlim_q50 = NA, pBlim_q05 = NA,
                            pBpa_q95 = NA,  pBpa_q50 = NA,  pBpa_q05 = NA,
                            pBtarget_q95 = NA, pBtarget_q50 = NA, pBtarget_q05 = NA,
                            pFlim_q95 = NA, pFlim_q50 = NA, pFlim_q05 = NA,
                            pFpa_q95 = NA,  pFpa_q50 = NA,  pFpa_q05 = NA,
                            pFtarget_q95 = NA, pFtarget_q50 = NA, pFtarget_q05 = NA)
      
      resRP <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season, .data$stock) %>%  
        dplyr::summarise(pBlim_q50 = sum(lowerBlim)/dplyr::n(),
                         pBpa_q50  = sum(lowerBpa)/dplyr::n(),
                         pBtarget_q50 = sum(lowerBtarget)/dplyr::n(),
                         pFlim_q50 = sum(higherFlim)/dplyr::n(),
                         pFpa_q50  = sum(higherFpa)/dplyr::n(),
                         pFtarget_q50 = sum(higherFtarget)/dplyr::n())
    }
    else{
      
      res <- obj %>% dplyr::group_by(.data$scenario, .data$stock, .data$year) %>%  
        summarise_at(c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards', 'catch.iyv', 'land.iyv', 'disc.iyv'),
                     .funs =  p_funs)
      
      # RP indicator
      res <- res %>% mutate(pBlim_q95 = NA, pBlim_q50 = NA, pBlim_q05 = NA,
                            pBpa_q95 = NA,  pBpa_q50 = NA,  pBpa_q05 = NA,
                            pBtarget_q95 = NA, pBtarget_q50 = NA, pBtarget_q05 = NA,
                            pFlim_q95 = NA, pFlim_q50 = NA, pFlim_q05 = NA,
                            pFpa_q95 = NA,  pFpa_q50 = NA,  pFpa_q05 = NA,
                            pFtarget_q95 = NA, pFtarget_q50 = NA, pFtarget_q05 = NA)
      
      resRP <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$stock) %>%  
        dplyr::summarise(pBlim_q50 = sum(lowerBlim)/dplyr::n(),
                         pBpa_q50  = sum(lowerBpa)/dplyr::n(),
                         pBtarget_q50 = sum(lowerBtarget)/dplyr::n(),
                         pFlim_q50 = sum(higherFlim)/dplyr::n(),
                         pFpa_q50  = sum(higherFpa)/dplyr::n(),
                         pFtarget_q50 = sum(higherFtarget)/dplyr::n())
      
    }
    res$pBlim_q50    <- resRP$pBlim_q50
    res$pBpa_q50     <- resRP$pBpa_q50
    res$pBtarget_q50 <- resRP$pBtarget_q50
    res$pFlim_q50    <- resRP$pFlim_q50
    res$pFpa_q50     <- resRP$pFpa_q50
    res$pFtarget_q50 <- resRP$pFtarget_q50
    
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

fltSum <- function (obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = FALSE, InterestRate = 0.03,scenario = 'bc')
{
  
  # For avoiding warnings in R CMD CHECK
  grossValue <- costs <- salaries <- grossSurplus <- nVessels <- fep <- NULL
  
  fleets <- obj$fleets
  covars <- obj$covars
  
  #  fleets <- lapply(fleets, setUnitsNA)
  
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')
  
  if (flnms[1] == "all") flnms <- names(fleets)
  
  Dim <- dim(fleets[[1]]@effort[,years , ])[c(2, 4, 6)]
  Dimnm <- dimnames(fleets[[1]]@effort[,years , ])
  
  FLQ0.dimfleets <- FLQuant(0, dimnames = c(fleet = list(names(fleets)), dimnames(fleets[[1]]@effort[, , ])[2:6]))
  
  if(is.null(covars$Depreciation))covars$Depreciation <- FLQ0.dimfleets
  if(is.null(covars$Salaries))    covars$Salaries <- FLQ0.dimfleets   
  if(is.null(covars$MaxDays))     covars$MaxDays <- FLQuant(365/dim(fleets[[1]]@effort[, years, ])[4], dimnames = dimnames(FLQ0.dimfleets))
  if(is.null(covars$NumbVessels)) covars$NumbVessels <- FLQ0.dimfleets
  if(is.null(covars$CapitalCost)) covars$CapitalCost <- FLQ0.dimfleets
  
  if(byyear == F){
    
    res <- NULL
    res.fl <- NULL
    
    year = rep(years, prod(Dim[2:3])) 
    season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3])
    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), 1)   
    
    
    for (f in flnms) {
      
 #     f <- names(fleets) #loop
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      fleet <- rep(f, each = prod(Dim))
      
      temp.catch <- lapply(catchNames(fl), function(x) quantSums(unitSums(catchWStock.f(fl, x))))
      temp.landings <- lapply(catchNames(fl), function(x) quantSums(unitSums(landWStock.f(fl, x))))
      temp.discards <- lapply(catchNames(fl), function(x) quantSums(unitSums(discWStock.f(fl, x))))
      

      # dga: I don't know why the code below failed during the mixfish wg and I've to replace 
      # the code by the code below.
      # res.fl <- bind_cols(year=year, season=season,fleet=fleet, iter=iter,
      #                     catch=c(Reduce('+',temp.catch)[,years]),
      #                     landings=c(Reduce('+',temp.landings)[,years]),
      #                     discards=c(Reduce('+',temp.discards)[,years]),
      #                     capacity=c(fl@capacity[,years, ]),
      #                     effort=c(fl@effort[,years, ]),
      #                     fcosts=c(totfcost_flbeia(fl, covars, f)[,years, ]),
      #                     vcosts=c(totvcost_flbeia(fl)[,years, ]),
      #                     costs=c(costs_flbeia(fl, covars, f)[,years, ]),
      #                     fcosts=c(totfcost_flbeia(fl, covars, f)[,years, ]),
      #                     vcosts=c(totvcost_flbeia(fl)[,years, ]),
      #                     costs=c(costs_flbeia(fl, covars, f)[,years, ]),
      #                     grossValue=c(revenue_flbeia(fl)[,years, ]),
      #                     nVessels = c(covars[['NumbVessels']][f,years])) %>% 
      
      res.fl <- tibble(data.frame(year=year, season=season,fleet=fleet, iter=iter,
                            catch=c(Reduce('+',temp.catch)[,years]),
                            landings=c(Reduce('+',temp.landings)[,years]),
                            discards=c(Reduce('+',temp.discards)[,years]),
                            capacity=c(fl@capacity[,years, ]),
                            effort=c(fl@effort[,years, ]),
                            fcosts=c(totfcost_flbeia(fl, covars, f)[,years, ]),
                            vcosts=c(totvcost_flbeia(fl)[,years, ]),
                            costs=c(costs_flbeia(fl, covars, f)[,years, ]),
                            grossValue=c(revenue_flbeia(fl)[,years, ]),
                            nVessels = c(covars[['NumbVessels']][f,years]))) %>% 
        mutate(discRat=discards/catch,
               grossSurplus=grossValue-costs,
               price=grossValue/landings,
               salaries=c(fl@crewshare[,years,])*grossValue+c(covars[['Salaries']][f,years]),
               gva=grossValue-costs+salaries,
               profitability=grossSurplus/grossValue,
               fep=grossSurplus- c(covars[['Depreciation']][f,years])*nVessels ,
               netProfit=fep- c(covars[['CapitalCost']][f,years])*InterestRate*nVessels)
      
      #quotaUptake depends on the number of seasons
      
      temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(catchWStock.f(fl, x))))
      temp <- Reduce('+',temp)[,years]
      totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,years]*obj$advice$TAC[x,years]))
      if (dim(temp)[4] >1) {
        res.fl <- res.fl %>% mutate(quotaUpt = c(sweep(temp, c(1:3,5:6),totTAC/dim(temp)[4], "/")))
      } else {
        res.fl <- res.fl %>% mutate(quotaUpt =  c(temp/totTAC))
      }
      res <- bind_rows(res,res.fl)
    }
  }else{
    
    res <- NULL
    res.fl <- NULL
    
    year = rep(years, Dim[3]) 
    iter = rep(rep(1:Dim[3], each = Dim[1]), 1)   
    
    
    for (f in flnms) {
      
 #     f <- names(fleets) #loop
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      fleet <- rep(f, each = prod(Dim[-2]))
      
      temp.catch <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
      temp.landings <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(landWStock.f(fl, x)))))
      temp.discards <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(discWStock.f(fl, x)))))
      

      # dga: I don't know why the code below failed during the mixfish wg and I've to replace 
      # the code by the code below.
      # res.fl <- bind_cols(year=year,fleet=fleet, iter=iter,
      #                     catch=c(Reduce('+',temp.catch)[,years]),
      #                     landings=c(Reduce('+',temp.landings)[,years]),
      #                     discards=c(Reduce('+',temp.discards)[,years]),
      #                     capacity=c(seasonSums(fl@capacity[,years, ])),
      #                     effort=c(seasonSums(fl@effort[,years, ])),
      #                     fcosts=c(seasonSums(totfcost_flbeia(fl, covars, f)[,years, ])),
      #                     vcosts=c(seasonSums(totvcost_flbeia(fl)[,years, ])),
      #                     costs=c(seasonSums(costs_flbeia(fl, covars, f)[,years, ])),
      #                     fcosts=c(seasonSums(totfcost_flbeia(fl, covars, f)[,years, ])),
      #                     vcosts=c(seasonSums(totvcost_flbeia(fl)[,years, ])),
      #                     costs=c(seasonSums(costs_flbeia(fl, covars, f)[,years, ])),
      #                     grossValue=c(seasonSums(revenue_flbeia(fl)[,years, ])),
      #                     nVessels =c(seasonMeans(covars[['NumbVessels']][f,years]))) %>%
        res.fl <- tibble(data.frame(year=year,fleet=fleet, iter=iter,
                            catch=c(Reduce('+',temp.catch)[,years]),
                            landings=c(Reduce('+',temp.landings)[,years]),
                            discards=c(Reduce('+',temp.discards)[,years]),
                            capacity=c(seasonSums(fl@capacity[,years, ])),
                            effort=c(seasonSums(fl@effort[,years, ])),
                            fcosts=c(seasonSums(totfcost_flbeia(fl, covars, f)[,years, ])),
                            vcosts=c(seasonSums(totvcost_flbeia(fl)[,years, ])),
                            costs=c(seasonSums(costs_flbeia(fl, covars, f)[,years, ])),
                            grossValue=c(seasonSums(revenue_flbeia(fl)[,years, ])),
                            nVessels =c(seasonMeans(covars[['NumbVessels']][f,years])))) 

               
        res.fl <- res.fl %>%  mutate(discRat = discards/catch,
               grossSurplus = grossValue - costs,
               price = grossValue/landings,
               salaries = c(seasonSums(fl@crewshare[,years,])) * grossValue + c(seasonSums(covars[['Salaries']][f,years])),
               gva = grossValue - costs + salaries,
               profitability = grossSurplus/grossValue,
               fep = grossSurplus - c(seasonSums(covars[['Depreciation']][f,years])) * nVessels ,
               netProfit = fep - c(seasonSums(covars[['CapitalCost']][f,years])) * InterestRate * nVessels)
      
      #quotaUptake depends on the number of seasons
      
      temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
      temp <- Reduce('+',temp)[,years]
      totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,years]*obj$advice$TAC[x,years]))
      res.fl <- res.fl %>% mutate(quotaUpt=c(temp/totTAC))
      res<- bind_rows(res,res.fl)
    }
  }
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "year"        "fleet"    "iter"       "catch":"quotaUpt"
  #  numeric      factor     numeric      numeric
  
  res <- res %>% mutate( year = as.numeric(year),
                         iter = as.numeric(iter))
  
  if(long == TRUE){ # transform res into long format
    ind <- if_else(byyear == TRUE , 3,4)
    indicator.nms <- names(res)[-c(1:ind)]
    res <- res %>% gather(key='indicator',value='value',indicator.nms)
    
  }
  
  res <- res  %>% mutate(scenario=scenario) %>% select(scenario, everything())
  return(res)
}

#' @rdname bioSum
#' @aliases fltSumQ

fltSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  p_names <- paste("q",ifelse(nchar(substr(prob,3, nchar(prob)))==1, 
                              paste(substr(prob,3, nchar(prob)), 0, sep = ""), 
                              substr(prob,3, nchar(prob))), sep = "")
  
  
  if(dim(obj)[2] < 10){ # the object is in long format
    if('season' %in% names(obj)){
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year,.data$season,.data$fleet,.data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
    }else{
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year,.data$fleet,.data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')}
    
  }else{
    p_funs <- purrr::map(prob, ~purrr::partial(quantile, prob = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    if('season' %in% names(obj)){
      sum.nms <- names(obj)[-c(1:5)]
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year,.data$season,.data$fleet) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)
    }
    else{
      sum.nms <- names(obj)[-c(1:4)]
      res <- obj %>% dplyr::group_by(.data$scenario, .data$year,.data$fleet) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)      
    }}
  
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
# @inheritParams FLBEIA
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
#        landings, discards, catch, discRat, price, tacshare, quota, quotaUpt, choke] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases fltStkSum
fltStkSum <- function(obj, flnms = names(obj$fleets), 
  stknms = catchNames(obj$fleets), 
  years = dimnames(obj$biols[[1]]@n)[[2]], 
  byyear = TRUE, long = FALSE, scenario = 'bc', 
  verbose = TRUE){
  
  fleets <- obj$fleets
  advice <- obj$advice
  
  fleets.ctrl <- obj$fleets.ctrl
  
  #  fleets <- lapply(fleets, setUnitsNA)

  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend
          removing the units using the setUnitsNA function')
  
  if(flnms[1] == 'all') flnms <- names(fleets)
  if(stknms[1] == 'all') stknms <- catchNames(fleets)
  
  Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
  Dimnm <- dimnames(fleets[[1]]@effort[,years,])
  
  resfl <- vector("list", length(flnms))
  names(resfl) <- flnms
  if(byyear == FALSE){ 
    
    year = rep(years, prod(Dim[2:3])) 
    season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3])
    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), 1)   
    
    
    for(f in flnms){
      
      fl   <- fleets[[f]]
      stfl <- catchNames(fl)        
      sts  <- stknms[stknms %in% stfl]
      
      resflst <- vector("list", length(sts))
      names(resflst) <- sts
      for(st in sts){
 
        # fl <- fleets[[f]]
        fleet = rep(f, each = prod(Dim))
        stock = rep(st, each = prod(Dim))
        
        tacshare <- sweep(fleets.ctrl$seasonal.share[[st]][f,], c(1:3,5:6), advice$quota.share[[st]][f,], '*')
        quota    <- sweep(tacshare, c(1:3,5:6), advice$TAC[st,], '*')
        
        # # checks
        # seasonSums(tacshare) == advice$quota.share[[st]][f,]
        # seasonSums(quota)    == advice$TAC[st,] * advice$quota.share[[st]][f,]
        
        res.fl.st <- bind_cols(year=year, season=season,fleet=fleet, stock=stock,iter=iter,
                               catch= c(apply(catchWStock.f(fl, st),c(2,4,6), sum)[,years]),
                               landings= c(apply(landWStock.f(fl, st),c(2,4,6), sum)[,years]),
                               discards= c(apply(discWStock.f(fl, st),c(2,4,6), sum)[,years])) %>% 
          mutate(discRat=discards/catch,
                 price=c(price_flbeia(fl, st)[,years]),
                 tacshare=c(tacshare[,years]),
                 quota=c(quota[,years]),
                 quotaUpt=catch/quota, choke = ifelse(round(catch/quota,2)>=1, TRUE, FALSE))
        resflst[[st]] <- res.fl.st
        if(verbose){print(paste("| fleet =", f, "|", "stock =", st, "|"))}
      }
      resfl[[f]] <- do.call("rbind", resflst)
    }
    res <- do.call("rbind", resfl)
  }  else{
    
    year = rep(years, Dim[3]) 
    iter = rep(rep(1:Dim[3], each = Dim[1]), 1)   
    
    for(f in flnms){
      
      fl   <- fleets[[f]]
      stfl <- catchNames(fl)        
      sts   <- stknms[stknms %in% stfl]
      fleet = rep(f, each = prod(Dim[-2]))
      
      resflst <- vector("list", length(sts))
      names(resflst) <- sts
      
      for(st in sts){
        
        stock = rep(st, each = prod(Dim[-2]))
        
        res.fl.st <- bind_cols(year=year,fleet=fleet, stock=stock, iter=iter,
                               catch= c(apply(catchWStock.f(fl, st),c(2,6), sum)[,years]),
                               landings= c(apply(landWStock.f(fl, st),c(2,6), sum)[,years]),
                               discards= c(apply(discWStock.f(fl, st),c(2,6), sum)[,years])) %>% 
          mutate(discRat=discards/catch,
                 price=c(seasonMeans(price_flbeia(fl, st)[,years]*quantSums(unitSums(landWStock.f(fl, st)[,years])))/landings),
                 tacshare=c(advice$quota.share[[st]][f,][,years]),
                 quota=c((advice$TAC[st,]*advice$quota.share[[st]][f,])[,years]),
                 quotaUpt=catch/quota, choke = ifelse(round(catch/quota,2)>=1, TRUE, FALSE))
        resflst[[st]] <- res.fl.st
        if(verbose){print(paste("| fleet =", f, "|", "stock =", st, "|"))}
      }
      resfl[[f]] <- do.call("rbind", resflst)
    }
    res <- do.call("rbind", resfl)
  }
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "year"     "fleet"    "stock"    "iter"     "catch":"quotaUpt"
  # numeric    character  character  numeric    numeric
  
  res <- res %>% mutate( year = as.numeric(year),
                         iter = as.numeric(iter))
  
  if(long == TRUE){ # transform res into long format
    ind <- if_else(byyear == TRUE , 4,5)
    indicator.nms <- names(res)[-c(1:ind)]
    res <- res %>% gather(key='indicator',value='value',indicator.nms)
  }
  
  res <- res  %>% mutate(scenario=scenario) %>% select(scenario, everything())
  
  return(res)
} 


# fltStkSumQ 
#' @rdname bioSum
#' @aliases fltStkSumQ


fltStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  
  
  p_names <- paste("q",ifelse(nchar(substr(prob,3, nchar(prob)))==1, 
                              paste(substr(prob,3, nchar(prob)), 0, sep = ""), 
                              substr(prob,3, nchar(prob))), sep = "")
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    # The choke indicator is treated differently, as it is a logical variable, then when doing summary we calculate the proportion
    # and not the quantiles.
    objCh <- obj %>% filter(indicator == 'choke')
    obj   <- obj %>% filter(indicator != 'choke')
    
    if('season' %in% names(obj)){
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season, .data$fleet, .data$stock, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
      # choke indicator
      resCh <- objCh %>% dplyr::group_by(.data$scenario,.data$year,.data$season, .data$fleet, .data$stock, .data$indicator) %>%
        dplyr::summarise(value = sum(value, na.rm=T)/dplyr::n()) 
      
      resCh <- bind_cols(resCh[,1:6],NA, resCh[,7], NA)
      names(resCh)[7:9] <- names(res)[7:9]
    }else{
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year,.data$fleet, .data$stock, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value))  %>% tidyr::spread(key='quantiles', value='value')
      # choke indicator
      resCh <- objCh %>% dplyr::group_by(.data$scenario,.data$year,.data$fleet, .data$stock, .data$indicator) %>%
        dplyr::summarise(value = sum(value, na.rm=T)/dplyr::n()) 
      
      resCh <- bind_cols(resCh[,1:5],NA, resCh[,6], NA)
      names(resCh)[6:8] <- names(res)[6:8]
    }
    res   <- bind_rows(res, resCh)
    
  }else{ # the object is in wide format
    p_funs <- purrr::map(prob, ~purrr::partial(quantile, prob = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    if('season' %in% names(obj)){
      sum.nms <- names(obj)[-c(1:6)]
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season,.data$fleet, .data$stock) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)
      
      # choke indicator
      res <- res %>% mutate(choke_q95 = NA, choke_q50 = NA, choke_q05 = NA)
      resCh <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season,.data$fleet, .data$stock) %>%  
        dplyr::summarise(choke_q50 = sum(choke, na.rm=T)/dplyr::n()) 
      res$choke_q50 <- resCh$choke_q50
    }
    else{
      sum.nms <- names(obj)[-c(1:5)]
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$fleet, .data$stock) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)   
      # choke indicator
      res <- res %>% mutate(choke_q95 = NA, choke_q50 = NA, choke_q05 = NA)
      resCh <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$fleet, .data$stock) %>%  
        dplyr::summarise(choke_q50 = sum(choke, na.rm=T)/dplyr::n()) 
      res$choke_q50 <- resCh$choke_q50
    }}
  
  
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
                     years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = FALSE, scenario = 'bc'){
  
  # for avoiding a warning in R CMD CHECK
  revst <- NULL
  
  
  fleets <- obj$fleets
  advice <- obj$advice
  
  # fleets <- lapply(fleets, setUnitsNA)
  
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')
  
  if(flnms[1] == 'all') flnms <- names(fleets)
  if(stknms[1] == 'all') stknms <- catchNames(fleets)
  
  
  Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
  Dimnm <- dimnames(fleets[[1]]@effort[,years,])
  
  res <- NULL
  res.fl.mt.stk <- NULL
  
  if(byyear == FALSE){
    
    year = rep(years, prod(Dim[2:3])) 
    season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3])
    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), 1)   
    
    for(f in flnms){
      
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      fleet = rep(f, each = prod(Dim))
      
      
      for(m in mts){
        
        mt <- fl@metiers[[m]]
        metier = rep(m, each = prod(Dim))
        
        stmt <- catchNames(mt)        
        sts  <- stknms[stknms %in% stmt]
        
        for(ss in sts){
          
          stock = rep(ss, each = prod(Dim))
          
          cc <- mt@catches[[ss]]
          
          res.fl.mt.ss <- bind_cols(year=year, season=season,fleet=fleet,metier=metier,stock=stock, iter=iter,
                                    landings=c(apply(cc@landings[,years,], c(2,4,6), sum,  na.rm=TRUE)),
                                    discards=c(apply(cc@discards[,years,], c(2,4,6), sum,  na.rm=TRUE)),
                                    revst =c(apply(cc@landings.n*cc@landings.wt*cc@price, c(2,4,6), sum, na.rm=TRUE)[,years,])) %>% 
            mutate(catch=landings+discards,
                   discRat=discards/catch,
                   price=revst/landings)
          
          res <- bind_rows(res, res.fl.mt.ss)           
          
        }
      }  
    }}
  else {
    year = rep(years, Dim[3]) 
    iter = rep(rep(1:Dim[3], each = Dim[1]), 1)   
    
    for(f in flnms){
      
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      fleet = rep(f, each = prod(Dim[-2]))
      
      for(m in mts){
        
        mt <- fl@metiers[[m]]
        metier = rep(m, each = prod(Dim[-2]))
        
        stmt <- catchNames(mt)        
        sts  <- stknms[stknms %in% stmt]
        
        
        for(ss in sts){
          
          stock = rep(ss, each = prod(Dim[-2]))
          
          cc <- mt@catches[[ss]]
          
          res.fl.mt.ss <- bind_cols(year=year,fleet=fleet,metier=metier,stock=stock, iter=iter,
                                    landings=c(apply(cc@landings[,years,], c(2,6), sum,  na.rm=TRUE)),
                                    discards=c(apply(cc@discards[,years,], c(2,6), sum,  na.rm=TRUE)),
                                    revst =c(apply(cc@landings.n*cc@landings.wt*cc@price, c(2,6), sum, na.rm=TRUE)[,years,])) %>% 
            mutate(catch=landings+discards,
                   discRat=discards/catch,
                   price=revst/landings)
          
          res <- bind_rows(res, res.fl.mt.ss)           
          
        }
      }  
    }
    
  }
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "year"     "fleet"    "metier"   "stock"    "iter"     "landings":"price"
  # numeric    character  character  character  numeric    numeric
  
  res <- res %>% mutate( year = as.numeric(year),
                         iter = as.numeric(iter))
  
  if(long == TRUE){ # transform res into long format
    ind <- if_else(byyear == TRUE , 5,6)
    indicator.nms <- names(res)[-c(1:ind)]
    res <- res %>% gather(key='indicator',value='value',indicator.nms)
  }
  
  res <- res  %>% mutate(scenario=scenario) %>% select(scenario, everything())
  
  return(res)
  
}

# mtStkSumQ 
#' @rdname bioSum
#' @aliases mtStkSumQ
mtStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  p_names <- paste("q",ifelse(nchar(substr(prob,3, nchar(prob)))==1, 
                              paste(substr(prob,3, nchar(prob)), 0, sep = ""), 
                              substr(prob,3, nchar(prob))), sep = "")
  
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if('season' %in% names(obj)){
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season, .data$fleet,.data$metier, .data$stock, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
    }else{
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year,.data$fleet,.data$metier, .data$stock, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value))  %>% tidyr::spread(key='quantiles', value='value')}
  }else{
    p_funs <- purrr::map(prob, ~purrr::partial(quantile, prob = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    if('season' %in% names(obj)){
      sum.nms <- names(obj)[-c(1:7)]
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season,.data$fleet,.data$metier, .data$stock) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)
    }
    else{
      sum.nms <- names(obj)[-c(1:6)]
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$fleet, .data$metier, .data$stock) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)      
    }}
  
  return(res)
}

#------------------------------------------------------------------------------#
# mtSum data.frame[scenario, year, season, fleet, metier, iter ||,|| 
#        effshare, effort, grossValue, vcost] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases mtSum

mtSum <- function(obj, flnms = names(obj$fleets),
                  years = dimnames(obj$biols[[1]]@n)[[2]], 
                  byyear = TRUE, long = FALSE, scenario = 'bc'){
  
  fleets <- obj$fleets
  
  #  fleets <- lapply(fleets, setUnitsNA)
  
  warning('Due to a problem with the units attribute in some off the slots, sometimes this function crashes. In case it fails, we recommend 
          removing the units using the setUnitsNA function')
  
  if (flnms[1] == "all") flnms <- names(fleets)
  
  
  Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
  Dimnm <- dimnames(fleets[[1]]@effort[,years,])
  
  res <- NULL
  res.fl.mt <- NULL
  
  
  if(byyear == FALSE){  
    
    year = rep(years, prod(Dim[2:3])) 
    season = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3])
    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), 1)   
    
    
    for(f in flnms){
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      fleet = rep(f, each = prod(Dim))
      
      
      for(m in mts){
        mt <- fl@metiers[[m]]
        metier = rep(m, each = prod(Dim))
        
        res.fl.mt <- bind_cols(year=year, season=season,fleet=fleet,metier=metier, iter=iter,
                               effshare=c(mt@effshare[,years,])) %>% 
          mutate(effort=c((fl@effort[,years,])*effshare),
                 vcost=c(mt@vcost[,years,])*effort,
                 grossValue=c(Reduce('+', lapply(mt@catches, 
                                                 function(x) unitSums(quantSums(x@landings.n*x@price))[,years]))))
        
        res <- bind_rows(res, res.fl.mt)           
      }
    }
  } else{
    
    
    year = rep(years, Dim[3]) 
    iter = rep(rep(1:Dim[3], each = Dim[1]), 1)   
    
    for(f in flnms){
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      fleet = rep(f, each = prod(Dim[-2]))
      
      
      for(m in mts){
        mt <- fl@metiers[[m]]
        metier = rep(m, each = prod(Dim[-2]))
        
        res.fl.mt <- bind_cols(year=year,  fleet=fleet,metier=metier,iter=iter,
                               effshare=c(seasonSums(mt@effshare[,years,]))) %>% 
          mutate(effort=c(seasonSums((fl@effort*mt@effshare)[,years,])),
                 vcost=c(seasonSums(mt@vcost[,years,]))*effort,
                 grossValue=c(Reduce('+', lapply(mt@catches,
                                                 function(x) seasonSums(unitSums(quantSums(x@landings.n*x@price)))[,years]))))          
        res <- bind_rows(res, res.fl.mt)           
      }
    }
  }
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "year"     "fleet"    "metier"   "iter"     "effshare":"grossValue"
  # numeric    character  character  numeric    numeric
  
  res <- res %>% mutate( year = as.numeric(year),
                         iter = as.numeric(iter))
  
  if(long == TRUE){ # transform res into long format
    ind <- if_else(byyear == TRUE , 4,5)
    indicator.nms <- names(res)[-c(1:ind)]
    res <- res %>% gather(key='indicator',value='value',indicator.nms)
  }
  
  res <- res  %>% mutate(scenario=scenario) %>% select(scenario, everything())
  
  return(res)
}

# mtStkSumQ 
#' @rdname bioSum
#' @aliases mtSumQ

mtSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  
  p_names <- paste("q",ifelse(nchar(substr(prob,3, nchar(prob)))==1, 
                              paste(substr(prob,3, nchar(prob)), 0, sep = ""), 
                              substr(prob,3, nchar(prob))), sep = "")
  
  
  if(dim(obj)[2] < 9){ # the object is in long format
    
    if('season' %in% names(obj)){
      res <- obj %>% dplyr::group_by(.data$scenario, .data$year, .data$season, .data$fleet, .data$metier, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
    }else{
      res <- obj %>% dplyr::group_by(.data$scenario, .data$year,.data$fleet, .data$metier, .data$indicator) %>%
        dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
        unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')}
  }else{
    p_funs <- purrr::map(prob, ~purrr::partial(quantile, prob = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    if('season' %in% names(obj)){
      sum.nms <- names(obj)[-c(1:6)]
      res <- obj %>% dplyr::group_by(.data$scenario,.data$year, .data$season,.data$fleet, .data$metier) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)
    }
    else{
      sum.nms <- names(obj)[-c(1:5)]
      res <- obj %>% dplyr::group_by(.data$scenario, .data$year, .data$fleet, .data$metier) %>%  
        summarise_at(c(sum.nms),.funs =  p_funs)      
    }}
  
  return(res)
}
#------------------------------------------------------------------------------#
# advSum :: data.frame[scenario, year, stock, iter ||,|| 
#        catch, discards, discRat, landings, quotaUpt, tac] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
#' @aliases advSum

advSum <- function(obj, stknms = 'all', years = dimnames(obj$biols[[1]]@n)$year, long = FALSE, scenario = 'bc'){
  
  # for avoiding a warning in R CMD CHECK
  tac <- NULL
  
  if(stknms == 'all') stknms <- names(obj$biols)  
  
  x1 <- Reduce(rbind, lapply(stknms, function(x)  cbind(stock = x, 
                                                               array2df(apply(catchWStock(obj$fleets, x), c(2,6), sum), label.x = 'catch')[,c('year', 'iter', 'catch')])))
  
  x2 <- Reduce(rbind, lapply(stknms, function(x)  cbind(stock = x, 
                                                               array2df(apply(landWStock(obj$fleets, x), c(2,6), sum), label.x = 'landings')[,c('year', 'iter', 'landings')])))
  
  x3 <- Reduce(rbind, lapply(stknms, function(x)  cbind(stock = x, 
                                                               array2df(apply(discWStock(obj$fleets, x), c(2,6), sum), label.x = 'discards')[,c('year', 'iter', 'discards')])))
  
  res <- as_tibble(cbind(x1,discards = x3[,4], landings = x2[,4]))
                
  x4 <- as_tibble(array2df(obj$advice$TAC, label.x = 'tac')[,c('stock', 'year', 'iter', 'tac')])
                
  res <- full_join(res, x4, by = c('stock', 'year', 'iter')) %>% mutate(scenario = scenario)
  
  # Wide format 
  res <- res %>%  dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$iter) %>% mutate(quotaUpt = catch/tac, discRat = discards/catch)
                
  # Reorder
  res <- res %>% select(stock, year, iter, scenario, !c('stock','year','iter','scenario'))
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "stock"     "year"      "iter"      "scenario"  "catch":"discRat"
  # character   numeric     numeric      character  numeric
  
  res <- res %>% mutate( year = as.numeric(as.character(year)),
                         iter = as.numeric(as.character(iter)))
  
  # Set consistent classess in summary outputs
  # lapply(res, class)
  # "stock"    "year"     "iter"     "catch":"discRat"
  # character   numeric   numeric    numeric
  
  res <- res %>% mutate( year = as.numeric(year),
                         iter = as.numeric(as.character(iter)))
  
  # reshaping this to the long format
   if(long == TRUE) res <- res %>% gather(key='indicator', value='value', .data$catch, .data$discards, .data$discRat, .data$landings, .data$quotaUpt, .data$tac)
                
 return(res)
}

#' @rdname bioSum
#' @aliases advSumQ
advSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  # for avoiding a warning in R CMD CHECK
  value <- NULL
  
  p_names <- paste("q",ifelse(nchar(substr(prob,3, nchar(prob)))==1, 
                              paste(substr(prob,3, nchar(prob)), 0, sep = ""), 
                              substr(prob,3, nchar(prob))), sep = "")
  
  if(dim(obj)[2] <= 7){ # the object is in long format
    
    res <- obj %>% dplyr::group_by(.data$scenario, .data$stock, .data$year, .data$indicator)  %>% 
      dplyr::summarise(quantiles = list(p_names), value=list(quantile(.data$value, prob=prob, na.rm = TRUE))) %>% 
      unnest(c(.data$quantiles,.data$value)) %>% tidyr::spread(key='quantiles', value='value')
    
  }
  else{
    
    p_funs <- purrr::map(prob, ~purrr::partial(quantile, prob = .x, na.rm = TRUE)) %>% 
      purrr::set_names(nm = p_names)
    
    res <- obj %>% dplyr::group_by(.data$scenario, .data$stock, .data$year) %>%  
      summarise_at(c("catch",    "discards", "discRat",  "landings", "quotaUpt", "tac"),
                   .funs =  p_funs)
    
  }
  
  return(res)
}

#----------------------------------------------------------------------
# riskSum(obj, stknms, Bpa, Blim, Prflim, flnms, years, scenario)
# Bpa = a named vector with the precautionary biomass per stock.
# Blim = a named vector with the limit biomass per stock.
# Prflim = a named vector with the limit profit per fleet.
#----------------------------------------------------------------------
#' @rdname bioSum
#' @aliases riskSum
riskSum <- function(obj, stknms = names(obj$biols), Bpa, Blim, Prflim, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], scenario = 'bc'){
  
  # For avoiding warnings in R CMD CHECK
  fleet <- grossSurplus <- refp <- year <- NULL
  
  
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
  
  bioS <- bioS %>% dplyr::group_by(.data$scenario, .data$year, .data$stock, .data$iter) %>% 
    mutate(Bpa = Bpa[stock], Blim = Blim[stock], risk.pa = as.numeric(ssb<Bpa), risk.lim = as.numeric(ssb<Blim))
  
  bioS.pa <- bioS %>% dplyr::group_by(.data$year, .data$stock, .data$scenario) %>% 
    dplyr::summarise(indicator="pBpa", value = sum(.data$risk.pa)/length(.data$risk.pa))
  
  bioS.lim <- bioS %>% dplyr::group_by(.data$year, .data$stock, .data$scenario) %>% 
    dplyr::summarise(indicator = "pBlim", value = sum(.data$risk.lim)/length(.data$risk.lim))
  
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
  
  flS <- flS %>% dplyr::group_by(.data$scenario, .data$year, .data$fleet, .data$iter) %>% 
    mutate(refp = Prflim[fleet], risk = as.numeric(grossSurplus < refp))
  
  outfl <- flS %>% dplyr::ungroup() %>%  dplyr::mutate(year = as.numeric(year)) %>% dplyr::group_by(.data$year, .data$fleet, .data$scenario) %>% 
    dplyr::summarise(indicator = "pPrflim", value = sum(.data$risk)/length(.data$risk)) %>% dplyr::rename(unit=fleet) 
  
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
  
  flS <- subset(flS, flS$year %in% c(years))
  
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
vesselSum <- function(obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = FALSE, scenario = 'bc'){
  
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
   res <- subset(flS, flS$indicator %in% c(ids, 'profitability'))
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
       res <- subset(flS, flS$indicator %in% c(ids, 'profitability'))
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
vesselStkSum <- function(obj, flnms = names(obj$fleets), stknms = catchNames(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = FALSE, scenario = 'bc'){
  
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
          sts <- catchNames(obj$fleets[[fl]])
          for(st in sts){
            flS[flS$fleet == fl & flS$stock == st,col] <- flS[flS$fleet == fl & flS$stock == st,col]/c(seasonMeans(covars[['NumbVessels']][fl,years]))
          }}}
      res <- flS
    }
    else{
      ids <- c("landings", "discards", "catch" ,   "price",  "tacshare",   "quota"   )
      for(id in ids){
        for(fl in flnms){
          sts <- catchNames(obj$fleets[[fl]])
          for(st in sts){
            flS[flS$indicator == id & flS$fleet == fl & flS$stock == st, 'value'] <- flS[flS$indicator == id & flS$fleet == fl & flS$stock == st, 'value']/c(seasonMeans(covars[['NumbVessels']][fl,years]))
      }}}
  
      res <- flS
    }}
  if(byyear == FALSE){
    if(long == FALSE){
      for(col in c(6:8,10:12)){
        for(fl in flnms){
          sts <- catchNames(obj$fleets[[fl]])
          for(st in sts){
            for(ss in dimnames(fleets[[fl]]@effort)[[4]]){
              flS[flS$fleet == fl & flS$stock == st & flS$season == ss,col] <- flS[flS$fleet == fl & flS$stock == st & flS$season == ss,col]/c((covars[['NumbVessels']][fl,years,,ss]))
          }}}}
      res <- flS
    }
    else{
      ids <- c("landings", "discards", "catch" ,   "price",  "tacshare",   "quota"   )
      for(id in ids){
        for(fl in flnms){
          sts <- catchNames(obj$fleets[[fl]])
          for(st in sts){
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


