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
# Dorleta GarcYYYa
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
#'       \item{F_flbeia}{ ithis function computes fishing mortiality.}
#'       \item{SSB_flbeia}{ this function computes spawning stock biomass by species.}
#'       \item{R_flbeia}{ this function computes recruitment by stock. If the stock is defined by age this function the recruiment is computed. ;
#'                        If the stock is follows a biomass dynamics, this function gives the growth.}
#'       \item{C_flbeia}{ this function computes catches by fleets and stock.} 
#'       \item{L_flbeia}{ this function computes landings by fleets and stock.}
#'       \item{D_flbeia}{ ithis function computes the discards by fleets and stock.}
#'      }     

#------------------------------------------------------------------------------#
# F_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
# @export

#' @rdname F_flbeia
F_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){
        # harvest: * if age structured calculate it from 'n'.
        #          * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
        na <- dim(obj$biols[[stk]]@n)[1]
        
        if(na == 1){
            # Catch:
            catch <- apply(catchStock(obj$fleets, stk),c(2,6), sum)[,years,drop = TRUE] # [ny,it]
            B     <- (obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,years,,1,drop= TRUE] # [ny, it] , 1st season biomass
            res[stk,,] <- (catch/B)
        }
        else{ 
            fbar_age <- ac(obj$biols[[stk]]@range[c('minfbar')]:obj$biols[[stk]]@range[c('maxfbar')])
            
            Dnms <- list(age = fbar_age, year = yrnms, iter = 1:it)
            aux  <- array(dim = c(length(fbar_age), ny,it), dimnames = Dnms)           
            
            n.  <- array(unitSums(obj$biols[[stk]]@n)[fbar_age,years,,1,drop=T], dim = c(length(fbar_age),ny,it), dimnames = Dnms)
            m.  <- array(seasonSums(unitMeans(obj$biols[[stk]]@m))[fbar_age,years,drop=T], dim = c(length(fbar_age),ny,it), dimnames = Dnms)
            c.  <- array(apply(catchStock(obj$fleets, stk),c(1:2,6), sum)[fbar_age,years, drop = TRUE], dim = c(length(fbar_age),ny,it), dimnames = Dnms)
        
            fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
        
            for(y in yrnms){
                for(a in fbar_age){
                    for(i in 1:it){
                        if(n.[a,y,i] == 0) aux[a,y,i] <- 0
                        else{
                           xx <- try(uniroot(fobj, lower = 0, upper = 1e6, n = n.[a,y,i], m=m.[a,y,i], c = c.[a,y,i])$root, silent = TRUE)
                           aux[a,y,i] <- ifelse(class(xx) == 'try-error', NA, xx)
                        }      
            }}}
           res[stk,,] <- apply(aux,2:3,mean) 
        }
    }
    return(res)
}


#------------------------------------------------------------------------------#
# SSB_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname F_flbeia
SSB_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years

    stknms <- names(obj$biols)

    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # SSB in spawning season
      # Before 2017: spawning season: first season with fraction of natural mortality before spawning < 1
      # Since 2017: SSB 1st January
      spwn.sson <- 1
      # si <- 0
      # while( (si-spwn.sson)!=0) { 
      #   si <- spwn.sson
      #   spwn.sson  <- ifelse( sum(spwn(obj$biols[[stk]])[ , , 1, spwn.sson, drop = T]<1,na.rm=T)==0, spwn.sson+1, spwn.sson)
      #   d  <- si-spwn.sson 
      # }
        res[stk,,] <- apply(unitSums(n(obj$biols[[stk]])*wt(obj$biols[[stk]])*fec(obj$biols[[stk]])*mat(obj$biols[[stk]]))[,years,,spwn.sson], c(2,6), sum, na.rm = TRUE)[drop=T]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# B_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname F_flbeia
B_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years

    stknms <- names(obj$biols)

    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(unitSums(obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,years,,1], c(2,6), sum, na.rm = TRUE)[drop=T]
    }
    return(res)
}

#------------------------------------------------------------------------------#
# R_flbeia(obj) :: res[stocks, years, it] 
# If age struc => recruitment.
# If biodyn    => growth.
#------------------------------------------------------------------------------#

#' @rdname F_flbeia
R_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years

    stknms <- names(obj$biols)

    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # 
        na <- dim(obj$biols[[stk]]@n)[1]
        # Recruitment season: first season with individuals at lower age class (Nage0>0)
        rec.sson <- 1
        si <- 0
        while( (si-rec.sson)!=0) { 
          si <- rec.sson
          rec.sson  <- ifelse( sum(obj$biols[[stk]]@n[1, , 1, rec.sson, drop = T]!=0,na.rm=T)==0, rec.sson+1, rec.sson)
          d  <- si-rec.sson 
        }
        if(na > 1){
            res[stk,,] <- obj$biols[[stk]]@n[1,years,1,rec.sson,drop=T]
            if(dim(obj$biols[[stk]]@n)[3]>1){
                for(ss in (rec.sson+1):dim(obj$biols[[stk]]@n)[3]) res[stk,,] <- res[stk,,] + obj$biols[[stk]]@n[1,years,ss,ss,drop=T]
            }
        }else{
            catch <- matrix(apply(catchStock(obj$fleets, stk),c(2,6), sum)[,years,drop = TRUE],ny,it) # [ny,it]
            B     <- matrix((obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,years,,1,drop= TRUE],ny,it) # [ny, it] , 1st season biomass
            res[stk,-ny,] <- B[-1,] - B[-ny,] + catch[-ny,]
            
        }
    }
    return(res)
}

#------------------------------------------------------------------------------#
# C_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname F_flbeia
C_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years

    stknms <- names(obj$biols)
    
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(catchWStock(obj$fleets, stk),c(2,6), sum)[,years,drop = TRUE] # [ny,it]
    }
    return(res)
}

#------------------------------------------------------------------------------#
# L_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
#' @rdname F_flbeia
L_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years

    stknms <- names(obj$biols)

    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(landWStock(obj$fleets, stk),c(2,6), sum)[,years,drop = TRUE] # [ny,it]
    }
    return(res)
}

#------------------------------------------------------------------------------#
# D_flbeia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#

#' @rdname F_flbeia
D_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    stknms <- names(obj$biols)

    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years)
    yrnms  <- years

    stknms <- names(obj$biols)

    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(discWStock(obj$fleets, stk),c(2,6), sum)[,years,drop = TRUE] # [ny,it]
    }
    return(res)
}

#' @rdname F_flbeia
summary_flbeia <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){

    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- length(years) # dim(obj$biols[[1]]@n)[2]
    yrnms  <- years # dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it, 7), dimnames = list(stock = stknms, year = yrnms, iter = 1:it, 
                                                      indicators = c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards')))
    
    res[,,,1] <- R_flbeia(obj,years)
    res[,,,2] <- SSB_flbeia(obj,years)
    res[,,,3] <- F_flbeia(obj,years)
    res[,,,4] <- B_flbeia(obj,years)
    res[,,,5] <- C_flbeia(obj,years)
    res[,,,6] <- L_flbeia(obj,years)
    res[,,,7] <- D_flbeia(obj,years)

    return(res)
    
}


#------------------------------------------------------------------------------#
# BIOsummary(obj) :: DATA.FRAME[stocks, years, it, indicators, value] 
#------------------------------------------------------------------------------#
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
#'              "fcosts", "gva", "income", "landings", "netProfit", "nVessels", "price", "profits",
#'              "quotaUpt", "salaries", "vcosts" and "profitability".
#'      \item{fltStkSum, fltStkSumQ:} Data frame with the indicators at fleet and stock level. The indicators are:
#'              "landings", "discards", "catch", "price",  "quotaUpt", "tacshare", "discRat" and  "quota".   
#'      \item{npv:} A data frame with the net present value per fleet over the selected range of years.
#'      \item{mtSum, mtSumQ:} Data frame with the indicators at fleet. The indicators are:
#'              "effshare", "effort", "income" and "vcost".   
#'      \item{mtStkSum, mtStkSumQ:} Data frame with the indicators at fleet and metier level. The indicators are:
#'              "catch",  "discards", "discRat", "landings" and "price".
#'      \item{riskSum:} A data frame with the risk indicators. The indicators are:
#'              "pBlim", "pBpa" and "pPrflim".
#'          
#'      
#'      
#'      \item{vesselSum, vesselSumQ:} Data frame with the indicators at vessel level. The indicators are:
#'               "catch", "costs", "discards", "discRat", "effort",       
#'              "fcosts", "gva", "income", "landings", "netProfit",  "price", "profits",
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
#' @return The data frames can be of wide or long format. In long format all the indicators are in the same column.
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
#' @inheritParams F_flbeia
#' @param flnms Names of the fleet for which the indicators will be calculated.
#' @param stknms Names of the stock for which the indicators will be calculated.
#' @param years the names of the years for which the indicators will be calculated. 
#' @param long logical. The data frame should be constructed using long or wide format? Default TRUE.
#' @param byyear logical. The indicators should be provided at season or year level? Default TRUE.
#' @param prob a numeric vector with the probabilities used to calculate the quantiles. 
#' @param scenario a character string with the name of the scenario corresponding with obj. Default bc.
#' @param Bpa named numeric vector with one element per stock in stknms. The precautionary approach stock spawning biomass used in riskSum function to calculate biological risk yearly.
#' @param Blim named numeric vector with one element per stock in stknms. The limit stock spawning biomass used in riskSum function to calculate biological risk yearly.
#' @param Prflim named numeric vector with one element per fleet in flnms. The limit profit level used in riskSum function to calculate economic risk yearly.

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
#' oneRes_bio    <- bioSum(oneRes, long = FALSE, years = ac(2016:2020))
#' oneRes_flt    <- fltSum(oneRes, long = FALSE, years = ac(2016:2020))
#' oneRes_fltStk <- fltStkSum(oneRes, long = FALSE, years = ac(2016:2020))
#' oneRes_mt     <- mtSum(oneRes, long = FALSE, years = ac(2016:2020))
#' oneRes_mtStk  <- mtStkSum(oneRes, long = FALSE, years = ac(2016:2020))
#' oneRes_adv    <- advSum(oneRes, long = FALSE, years = ac(2016:2020))
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
#' # Wide format with seasonal disaggregation. No seasonal disagregation available for bio and adv summaries.
#' 
#' oneRes_bio    <- bioSum(oneRes, long = FALSE) # Biol summary is only by year.
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
# 
#  # Long format and seasonal
#' oneRes_bio    <- bioSum(oneRes, long = TRUE) # Biol summary is only by year.
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
#' # Example OneIters: As one but with iterations.
#' #------------------------------------------------
#' s1_bio    <- bioSum(s1, scenario = 'with_iters')
#' s1_flt    <- fltSum(s1, scenario = 'with_iters')
#' s1_fltStk <- fltStkSum(s1, scenario = 'with_iters')
#' s1_mt     <- mtSum(s1, scenario = 'with_iters')
#' s1_mtStk  <- mtStkSum(s1, scenario = 'with_iters')
#' s1_adv    <- advSum(s1, scenario = 'with_iters')
#' 
#' s1_bioQ    <- bioSumQ(s1_bio)
#' s1_fltQ    <- fltSumQ(s1_flt)
#' s1_fltStkQ <- fltStkSumQ(s1_fltStk)
#' s1_mtQ     <- mtSumQ(s1_mt)
#' s1_mtStkQ  <- mtStkSumQ(s1_mtStk)
#' s1_advQ    <- advSumQ(s1_adv)
#' 
#' s1_bio    <- bioSum(s1, long = FALSE, years = ac(2016:2020))
#' s1_flt    <- fltSum(s1, long = FALSE, years = ac(2016:2020))
#' s1_fltStk <- fltStkSum(s1, long = FALSE, years = ac(2016:2020))
#' s1_mt     <- mtSum(s1, long = FALSE, years = ac(2016:2020))
#' s1_mtStk  <- mtStkSum(s1, long = FALSE, years = ac(2016:2020))
#' s1_adv    <- advSum(s1, long = FALSE, years = ac(2016:2020))
#' 
#' 
#' s1_bioQ    <- bioSumQ(s1_bio)
#' s1_fltQ    <- fltSumQ(s1_flt)
#' s1_fltStkQ <- fltStkSumQ(s1_fltStk)
#' s1_mtQ     <- mtSumQ(s1_mt)
#' s1_mtStkQ  <- mtStkSumQ(s1_mtStk)
#' s1_advQ    <- advSumQ(s1_adv)
#' 
#' 
#' s1_bio    <- bioSum(s1, long = FALSE) # Biol summary is only by year.
#' s1_flt    <- fltSum(s1, long = FALSE, byyear = FALSE)
#' s1_fltStk <- fltStkSum(s1, long = FALSE, byyear = FALSE)
#' s1_mt     <- mtSum(s1, long = FALSE, byyear = FALSE)
#' s1_mtStk  <- mtStkSum(s1, long = FALSE, byyear = FALSE)
#' s1_adv    <- advSum(s1, long = FALSE) # Advice summary is only by year.
#' 
#' s1_bioQ    <- bioSumQ(s1_bio)
#' s1_fltQ    <- fltSumQ(s1_flt)
#' s1_fltStkQ <- fltStkSumQ(s1_fltStk)
#' s1_mtQ     <- mtSumQ(s1_mt)
#' s1_mtStkQ  <- mtStkSumQ(s1_mtStk)
#' s1_advQ    <- advSumQ(s1_adv)
#' 
#' 
#' s1_bio    <- bioSum(s1, long = TRUE) # Biol summary is only by year.
#' s1_flt    <- fltSum(s1, long = TRUE, byyear = FALSE)
#' s1_fltStk <- fltStkSum(s1, long = TRUE, byyear = FALSE)
#' s1_mt     <- mtSum(s1, long = TRUE, byyear = FALSE)
#' s1_mtStk  <- mtStkSum(s1, long = TRUE, byyear = FALSE)
#' s1_adv    <- advSum(s1, long = TRUE) # Advice summary is only by year.
#' 
#' s1_bioQ    <- bioSumQ(s1_bio)
#' s1_fltQ    <- fltSumQ(s1_flt)
#' s1_fltStkQ <- fltStkSumQ(s1_fltStk)
#' s1_mtQ     <- mtSumQ(s1_mt)
#' s1_mtStkQ  <- mtStkSumQ(s1_mtStk)
#' s1_advQ    <- advSumQ(s1_adv)
#' 
#' s1_risk <- riskSum(s1, Bpa = c(stk1= 900), Blim = c(stk1 = 600), Prflim = c(flt1 = 0), scenario = 'alternative')
#' 
#' s1_npv  <- npv(s1, y0 = '2014')
#'
#' #------------------------------------------------
#' # Example Multi: Two stock, two fleet, four iters.
#' #------------------------------------------------
#' s2_bio    <- bioSum(s2)
#' s2_flt    <- fltSum(s2)
#' s2_fltStk <- fltStkSum(s2)
#' s2_mt     <- mtSum(s2)
#' s2_mtStk  <- mtStkSum(s2)
#' s2_adv    <- advSum(s2)
#' 
#' s2_bioQ    <- bioSumQ(s2_bio)
#' s2_fltQ    <- fltSumQ(s2_flt)
#' s2_fltStkQ <- fltStkSumQ(s2_fltStk)
#' s2_mtQ     <- mtSumQ(s2_mt)
#' s2_mtStkQ  <- mtStkSumQ(s2_mtStk)
#' s2_advQ    <- advSumQ(s2_adv)
#' 
#' s2_bio    <- bioSum(s2, long = FALSE, years = ac(2016:2020))
#' s2_flt    <- fltSum(s2, long = FALSE, years = ac(2016:2020))
#' s2_fltStk <- fltStkSum(s2, long = FALSE, years = ac(2016:2020))
#' s2_mt     <- mtSum(s2, long = FALSE, years = ac(2016:2020))
#' s2_mtStk  <- mtStkSum(s2, long = FALSE, years = ac(2016:2020))
#' s2_adv    <- advSum(s2, long = FALSE, years = ac(2016:2020))
#' 
#' 
#' s2_bioQ    <- bioSumQ(s2_bio)
#' s2_fltQ    <- fltSumQ(s2_flt)
#' s2_fltStkQ <- fltStkSumQ(s2_fltStk)
#' s2_mtQ     <- mtSumQ(s2_mt)
#' s2_mtStkQ  <- mtStkSumQ(s2_mtStk)
#' s2_advQ    <- advSumQ(s2_adv)
#' 
#' 
#' s2_bio    <- bioSum(s2, long = FALSE) # Biol summary is only by year.
#' s2_flt    <- fltSum(s2, long = FALSE, byyear = FALSE)
#' s2_fltStk <- fltStkSum(s2, long = FALSE, byyear = FALSE)
#' s2_mt     <- mtSum(s2, long = FALSE, byyear = FALSE)
#' s2_mtStk  <- mtStkSum(s2, long = FALSE, byyear = FALSE)
#' s2_adv    <- advSum(s2, long = FALSE) # Advice summary is only by year.
#' 
#' s2_bioQ    <- bioSumQ(s2_bio)
#' s2_fltQ    <- fltSumQ(s2_flt)
#' s2_fltStkQ <- fltStkSumQ(s2_fltStk)
#' s2_mtQ     <- mtSumQ(s2_mt)
#' s2_mtStkQ  <- mtStkSumQ(s2_mtStk)
#' s2_advQ    <- advSumQ(s2_adv)
#' 
#' 
#' s2_bio    <- bioSum(s2, long = TRUE) # Biol summary is only by year.
#' s2_flt    <- fltSum(s2, long = TRUE, byyear = FALSE)
#' s2_fltStk <- fltStkSum(s2, long = TRUE, byyear = FALSE)
#' s2_mt     <- mtSum(s2, long = TRUE, byyear = FALSE)
#' s2_mtStk  <- mtStkSum(s2, long = TRUE, byyear = FALSE)
#' s2_adv    <- advSum(s2, long = TRUE) # Advice summary is only by year.
#' 
#' s2_bioQ    <- bioSumQ(s2_bio)
#' s2_fltQ    <- fltSumQ(s2_flt)
#' s2_fltStkQ <- fltStkSumQ(s2_fltStk)
#' s2_mtQ     <- mtSumQ(s2_mt)
#' s2_mtStkQ  <- mtStkSumQ(s2_mtStk)
#' s2_advQ    <- advSumQ(s2_adv)
#' 
#' s2_npv  <- npv(s2, y0 = '2014')
#' risk_s2 <- riskSum(s2, Bpa = c(stk1= 135000, stk2 = 124000), Blim = c(stk1= 96000, stk2 = 89000), Prflim = c(flt1 = 0, flt2 = 0), scenario = 'alternative')
#'
#' }
bioSum <- function(obj, stknms = 'all', years = dimnames(obj$biols[[1]]@n)$year, long = TRUE, scenario = 'bc'){
    xx <- summary_flbeia(obj, years)
     
    dnms <- dimnames(xx)
    
    if(stknms[1] == 'all') stknms <- dnms[[1]]
    
    xx <- xx[stknms,,,,drop=F]
    dnms <- dimnames(xx)
    
    if(long == TRUE){
      
      df <- expand.grid(iter = dnms[[3]], year = dnms[[2]],indicator = c(dnms[[4]], 'land.iyv', 'disc.iyv', 'catch.iyv'),  stock = stknms)[,4:1]

      df$stock     <- as.character(df$stock)
      df$year      <- as.numeric(as.character(df$year))
      df$indicator <- as.character(df$indicator)
      df$iter      <- as.numeric(df$iter)

      df <- cbind(df, value = NA)
    
      for(st in dnms[[1]]){
            for(ind in dnms[[4]]){
               df[df$stock == st & df$indicator == ind,'value'] <- c(t(xx[st,,,ind]))
            }
            for(i in 1:length(dnms[[3]])){
              df[df$stock == st & df$indicator == 'land.iyv' & df$iter == i,'value'][-1] <- df[df$stock == st & df$indicator == 'landings' & df$iter == i,'value'][-1]/df[df$stock == st & df$indicator == 'landings' & df$iter == i,'value'][-length(years)]
              df[df$stock == st & df$indicator == 'disc.iyv' & df$iter == i,'value'][-1] <- df[df$stock == st & df$indicator == 'discards' & df$iter == i,'value'][-1]/df[df$stock == st & df$indicator == 'discards' & df$iter == i,'value'][-length(years)]
              df[df$stock == st & df$indicator == 'catch.iyv' & df$iter == i,'value'][-1] <- df[df$stock == st & df$indicator == 'catch' & df$iter == i,'value'][-1]/df[df$stock == st & df$indicator == 'catch' & df$iter == i,'value'][-length(years)]
              }
      }
      df <- df[,c('year','stock', 'iter', 'indicator','value')]
    }
    else{ # long = FALSE
      df <- expand.grid(iter = dnms[[3]], year = dnms[[2]],  stock = dnms[[1]])[,3:1]
      
      df$stock     <- as.character(df$stock)
      df$year      <- as.numeric(as.character(df$year))
      df$iter      <- as.numeric(df$iter)
      
      df <- cbind(df, biomass = NA, catch = NA, catch.iyv = NA,  discards = NA, disc.iyv = NA, f = NA, landings = NA, land.iyv = NA, rec = NA, ssb = NA)
      
      for(st in dnms[[1]]){
        for(ind in dnms[[4]]){
          df[df$stock == st,'biomass']  <- c(t(xx[st,,,'biomass'])) 
          df[df$stock == st,'catch']    <- c(t(xx[st,,,'catch']))  
          df[df$stock == st,'discards'] <- c(t(xx[st,,,'discards']))
          df[df$stock == st,'f']        <- c(t(xx[st,,,'f']))
          df[df$stock == st,'landings'] <- c(t(xx[st,,,'landings']))
          df[df$stock == st,'rec']      <- c(t(xx[st,,,'rec']))
          df[df$stock == st,'ssb']      <- c(t(xx[st,,,'ssb']))
          df[df$stock == st,'catch.iyv']     <- c(rep(NA, dim(xx)[3]), t(xx[st,-1,,'catch']/xx[st,-dim(xx)[2],,'catch']))
          df[df$stock == st,'land.iyv']      <- c(rep(NA, dim(xx)[3]), t(xx[st,-1,,'landings']/xx[st,-dim(xx)[2],,'landings']))
          df[df$stock == st,'disc.iyv']      <- c(rep(NA, dim(xx)[3]), t(xx[st,-1,,'discards']/xx[st,-dim(xx)[2],,'discards']))
          df <- df[,c('year','stock', 'iter', 'biomass', 'catch', 'discards', 
                          'f', 'landings', 'rec', 'ssb', 'catch.iyv', 'land.iyv', 'disc.iyv')]
          }
      }
    }
    
    df <- cbind(scenario = scenario, df)
    return(df)
}

#' @rdname bioSum
bioSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){

  if(dim(obj)[2] <= 6){ # the object is in long format
    res <- aggregate(value ~ stock + indicator + year + scenario, obj, quantile, prob = prob)
    res <- cbind(res[,1:4], data.frame(res[,5],stringsAsFactors = FALSE))
    
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
           
    names(res)[5:(5+length(prob)-1)] <- nms
  }
  else{
    res <- aggregate(list(biomass = obj$biomass,catch = obj$catch,catch.iyv = obj$catch.iyv, 
                          discards = obj$discards, disc.iyv = obj$disc.iyv,
                          f = obj$f,landings = obj$landings, land.iyv = obj$land.iyv,
                          rec = obj$rec, ssb = obj$ssb), 
                       list(stock = obj$stock, year = obj$year, scenario = obj$scenario), quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]), data.frame(res[,5]), data.frame(res[,6]),
                           data.frame(res[,7]), data.frame(res[,8]), data.frame(res[,9]),
                           data.frame(res[,10]), data.frame(res[,11]), data.frame(res[,12]),
                           data.frame(res[,13]))
                                                                                 
    nmsp <- ifelse(substr(prob,3, nchar(prob))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""),substr(prob,3, nchar(prob)))
    
    nms_bio  <- paste('biomass_q',nmsp, sep = "")
    nms_cat  <- paste('catch_q',nmsp, sep = "")
    nms_disc <- paste('discards_q',nmsp, sep = "")
    nms_land <- paste('landings_q',nmsp, sep = "")
    nms_f    <- paste('f_q',nmsp, sep = "")
    nms_rec  <- paste('rec_q',nmsp, sep = "")
    nms_ssb  <- paste('ssb_q',nmsp, sep = "")
    nms_iyvcat  <- paste('catch.iyv_q',nmsp, sep = "")
    nms_iyvdisc <- paste('disc.iyv_q',nmsp, sep = "")
    nms_iyvland <- paste('land.iyv_q',nmsp, sep = "")
    
    names(res)[-c(1:3)] <- c(nms_bio, nms_cat, nms_iyvcat, nms_disc, nms_iyvdisc, nms_f, nms_land, nms_iyvland,nms_rec, nms_ssb)
    
    
  }
  
  return(res)
  }
  
  
#------------------------------------------------------------------------------#
# fltSum data.frame[year, season, stock, fleet, iter, ||,|| 
#        profits, capacity, costs, discards, effort, landings] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
fltSum <- function (obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = TRUE, scenario = 'bc')
{
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
        income = numeric(n), 
        landings = numeric(n),
        netProfit = numeric(n), 
        nVessels = numeric(n), 
        price = numeric(n), 
        profits = numeric(n),
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
       
        res[k:(k + prod(Dim) - 1), "income"] <- c(revenue_flbeia(fl)[,years, ]) 
        
        res[k:(k + prod(Dim) - 1), "profits"] <- c(revenue_flbeia(fl)[,years, ]) - res[k:(k + prod(Dim) - 1), "costs"]
        
        res[k:(k + prod(Dim) - 1), "price"] <- res[k:(k + prod(Dim) - 1), "income"] / res[k:(k + prod(Dim) - 1), "landings"]
         
        res[k:(k + prod(Dim) - 1), "salaries"] <- c(fl@crewshare[,years,]*revenue_flbeia(fl)[,years, ] + covars[['Salaries']][f,years])
        
        res[k:(k + prod(Dim) - 1), "gva"] <- res[k:(k + prod(Dim) - 1), "income"] - res[k:(k + prod(Dim) - 1), "costs"] + res[k:(k + prod(Dim) - 1), "salaries"]
        
        res[k:(k + prod(Dim) - 1), "profitability"] <- res[k:(k + prod(Dim) - 1), "profits"]/res[k:(k + prod(Dim) - 1), "income"]
 
        res[k:(k + prod(Dim) - 1), "nVessels"]  <- c(covars[['NumbVessels']][f,years])
          
        res[k:(k + prod(Dim) - 1), "netProfit"] <- c(revenue_flbeia(fl)[,years, ] -  costs_flbeia(fl, covars, f)[,years, ] - covars[['Depreciation']][f,years]*covars[['NumbVessels']][f,years])

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
                        income = numeric(n), 
                        landings = numeric(n),
                        netProfit = numeric(n), 
                        nVessels = numeric(n), 
                        price = numeric(n), 
                        profits = numeric(n),
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
        
        res[k:(k + prod(Dim[-2]) - 1), "income"] <- c(seasonSums(revenue_flbeia(fl)[,years, ])) 
        
        res[k:(k + prod(Dim[-2]) - 1), "profits"] <- c(seasonSums(revenue_flbeia(fl)[,years, ])) - res[k:(k + prod(Dim[-2]) - 1), "costs"]
        
        res[k:(k + prod(Dim[-2]) - 1), "salaries"] <- c(seasonSums(fl@crewshare[,years,]*revenue_flbeia(fl)[,years, ] + covars[['Salaries']][f,years]))
        
        res[k:(k + prod(Dim[-2]) - 1), "gva"] <- res[k:(k + prod(Dim[-2]) - 1), "income"] - res[k:(k + prod(Dim[-2]) - 1), "costs"] + res[k:(k + prod(Dim[-2]) - 1), "salaries"]
        
        res[k:(k + prod(Dim[-2]) - 1), "profitability"] <- res[k:(k + prod(Dim[-2]) - 1), "profits"]/res[k:(k + prod(Dim[-2]) - 1), "income"]
        
        res[k:(k + prod(Dim[-2]) - 1), "nVessels"]  <- c(seasonMeans(covars[['NumbVessels']][f, years, ]))
        
        res[k:(k + prod(Dim[-2]) - 1), "netProfit"] <- c(seasonSums(revenue_flbeia(fl)[,years, ] -  costs_flbeia(fl, covars, f)[,years, ] - covars[['Depreciation']][f,years]*covars[['NumbVessels']][f,years]))
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
        temp <- Reduce('+',temp)[, years, ]
        totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,years]*obj$advice$TAC[x,years]))
        res[k:(k + prod(Dim[-2]) - 1), "quotaUpt"] <- c(totTAC/temp)
        
        k <- k + prod(Dim[-2])
    }}
    
    if(long == TRUE){ # transform res into long format
      r1 <- ifelse(byyear == TRUE, 4,5)
      r2 <- ifelse(byyear == TRUE, 21,22)
      
      names(res)[r1:r2] <- paste('indicator',names(res)[r1:r2], sep = "_")
      res <- reshape(res, direction = 'long', varying = r1:r2, sep = "_")[,1:(r1+1)]
      rownames(res) <- 1:dim(res)[1]
      names(res)[r1:(r1+1)] <- c('indicator', 'value') 
    }
  
   res <- cbind(scenario = scenario, res)
    return(res)
}
    
#' @rdname bioSum
fltSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
      res <- cbind(res[,1:4], data.frame(res[,5]))
    
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
      names(res)[5:(5+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + indicator + year + scenario + season, obj, quantile, prob = prob, na.rm=T)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[6:(6+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(capacity = obj$capacity,      catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                          discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                          income  = obj$income,         landings  = obj$landings,  netProfit  = obj$netProfit, nVessels  = obj$nVessels,     
                          price  = obj$price,           profits  = obj$profits,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                          vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                          list(fleet = obj$fleet, year = obj$year, scenario = obj$scenario), 
                          quantile, prob = prob, na.rm=T)
    
      res <- cbind(res[,1:3], 
                 data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                 data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                 data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                 data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]),
                 data.frame(res[,20]), data.frame(res[,21]))
                 
      nms1  <- paste('capacity_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('costs_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms8  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms9 <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms10 <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms11 <- paste('netprofit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms12 <- paste('nVessels_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms13 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms14 <- paste('profit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms15 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms16 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms17 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms18 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
   
      names(res)[-c(1:3)] <- unlist(mget(paste('nms', 1:18, sep="")))
    }
    else{
        res <- aggregate(list(capacity = obj$capacity,      catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                              discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                              income  = obj$income,         landings  = obj$landings,  netProfit  = obj$netProfit, nVessels  = obj$nVessels,     
                              price  = obj$price,           profits  = obj$profits,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                              vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                         list(fleet = obj$fleet, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.rm = TRUE)
        
        res <- cbind(res[,1:3], 
                     data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                     data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                     data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                     data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]),
                     data.frame(res[,20]), data.frame(res[,21]))
        
        nms1  <- paste('capacity_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms3  <- paste('costs_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms4  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms5  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms6  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms7  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms8  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms9 <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms10 <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms11 <- paste('netprofit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms12 <- paste('nVessels_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms13 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms14 <- paste('profit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms15 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms16 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms17 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms18 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        
        names(res)[-c(1:3)] <- unlist(mget(paste('nms', 1:18, sep="")))
      
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
            res <- res + FLQuant(apply(dat@landings.n*dat@landings.wt*dat@price, c(2,4,6),sum,na.rm=T), dim = D)
        }
    }
    return(res)               
}

#-------------------------------------------------------------------------------
# costs_flbeia(fleet, years)
#-------------------------------------------------------------------------------

#' @rdname revenue_flbeia
costs_flbeia <- function(fleet, covars, flnm = NULL){
    
    res <- totvcost_flbeia(fleet) + totfcost_flbeia(fleet, covars, flnm)
    
    return(res)               
}

#-------------------------------------------------------------------------------
# totvcost_flbeia(fleet, years)
#-------------------------------------------------------------------------------
#' @rdname revenue_flbeia
totvcost_flbeia <- function(fleet){
    
    mts <- names(fleet@metiers)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        res <- res + fleet@metiers[[mt]]@vcost*fleet@effort*fleet@metiers[[mt]]@effshare
    }
    Rev <- revenue_flbeia(fleet)*fleet@crewshare
    
    res <- res + Rev
    
    return(res)               
}

#-------------------------------------------------------------------------------
# totvcost_flbeia(fleet, years)
#-------------------------------------------------------------------------------
#' @rdname revenue_flbeia
totfcost_flbeia <- function(fleet, covars, flnm = NULL){
     if(is.null(flnm)) flnm <- 1
     return(fleet@fcost*covars[["NumbVessels"]][flnm, ])            
}


#------------------------------------------------------------------------------#
# catchFlSum data.frame[year, season, stock, fleet, iter, ||,|| 
#        landings, discards, price, tacshare] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
fltStkSum <- function(obj, flnms = names(obj$fleets), stknms = catchNames(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
    
  fleets <- obj$fleets
  advice <- obj$advice
  
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
                    price    = numeric(n),
                    quotaUpt = numeric(n),
                    tacshare = numeric(n),
                    stringsAsFactors = FALSE)
        
        k <- 1
        
        for(st in sts){
            
            dff[k:(prod(Dim) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,4,6), sum)[,years])    
            dff[k:(prod(Dim) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,4,6), sum)[,years]) 
            dff[k:(prod(Dim) + k-1),'catch']    <- dff[k:(prod(Dim) + k-1),'discards'] + dff[k:(prod(Dim) + k-1),'landings']
            dff[k:(prod(Dim) + k-1),'discRat']  <- dff[k:(prod(Dim) + k-1),'discards']/dff[k:(prod(Dim) + k-1),'catch']
            dff[k:(prod(Dim) + k-1),'price']    <- c(price_flbeia(fl, st)[,years])
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
        
        n <- prod(Dim)[-2]*length(sts)
        
        dff <- data.frame(year = rep(years, prod(Dim[3])*length(sts)), 
                          fleet = rep(f, n), 
                          stock = rep(sts, each = prod(Dim[-2])),
                          iter = rep(rep(1:Dim[3], each = prod(Dim[1])), length(sts)),  
                          landings = numeric(n), 
                          discards = numeric(n),
                          catch    = numeric(n),
                          price    = numeric(n),
                          quotaUpt = numeric(n),
                          tacshare = numeric(n),
                          stringsAsFactors = FALSE)
        
        k <- 1
        
        for(st in sts){
          
          dff[k:(prod(Dim[-2]) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,6), sum)[,years])    
          dff[k:(prod(Dim[-2]) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,6), sum)[,years]) 
          dff[k:(prod(Dim[-2]) + k-1),'catch']    <- dff[k:(prod(Dim[-2]) + k-1),'discards'] + dff[k:(prod(Dim[-2]) + k-1),'landings']
          dff[k:(prod(Dim[-2]) + k-1),'discRat']  <- dff[k:(prod(Dim[-2]) + k-1),'discards']/dff[k:(prod(Dim[-2]) + k-1),'catch']
          dff[k:(prod(Dim[-2]) + k-1),'price']    <- c(seasonMeans(price_flbeia(fl, st)[,years]*quantSums(unitSums(landWStock.f(fl, st)[,years])))/seasonSums(unitSums(quantSums(landWStock.f(fl, st)[,years]))))
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
fltStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[6:(6+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario + season, obj, quantile, prob = prob, na.rm=T)
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
                       quantile, prob = prob, na.rm=T)
      
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
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.rm = TRUE)
      
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



#-------------------------------------------------------------------------------
# price_flbeia(fleet, years)(mean price in a fleet)
#-------------------------------------------------------------------------------
#' @rdname revenue_flbeia
price_flbeia <- function(fleet, stock){

    mts <- names(fleet@metiers)
    
    totL <- apply(landWStock.f(fleet, stock), c(2,4,6), sum)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        m <- fleet@metiers[[mt]]
        if(!(stock %in% catchNames(m))) next
        dat <- m@catches[[stock]]
        res <- res + apply(dat@landings.n*dat@landings.wt*dat@price, c(2,4,6),sum,na.rm=T)
    }
    
    res <- res/totL
    
    return(res)                
}


#------------------------------------------------------------------------------#
# mtStkSum data.frame[year, season, stock, fleet, metier, iter ||,|| 
#        landings, discards, price] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
mtStkSum <- function(obj, flnms = names(obj$fleets), stknms = catchNames(obj$fleets), 
                     years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
    
    
  fleets <- obj$fleets
  advice <- obj$advice
  
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
                dfm[k:(k+prod(Dim)-1),'landings'] <- c(apply(cc@landings[,years,], c(2,4,6), sum, na.rm=T))
                dfm[k:(k+prod(Dim)-1),'discards'] <- c(apply(cc@discards[,years,], c(2,4,6), sum, na.rm=T))
                dfm[k:(k+prod(Dim)-1),'catch'] <- dfm[k:(k+prod(Dim)-1),'discards'] + dfm[k:(k+prod(Dim)-1),'landings']
                dfm[k:(k+prod(Dim)-1),'discRat'] <-  dfm[k:(k+prod(Dim)-1),'discards']/dfm[k:(k+prod(Dim)-1),'catch']
                revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,4,6), sum, na.rm=T)[,years,]
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
            dfm[k:(k+prod(Dim[-2])-1),'landings'] <- c(apply(cc@landings[,years,], c(2,6), sum, na.rm=T))
            dfm[k:(k+prod(Dim[-2])-1),'discards'] <- c(apply(cc@discards[,years,], c(2,6), sum, na.rm=T))
            dfm[k:(k+prod(Dim[-2])-1),'catch'] <- dfm[k:(k+prod(Dim[-2])-1),'discards'] + dfm[k:(k+prod(Dim[-2])-1),'landings']
            dfm[k:(k+prod(Dim[-2])-1),'discRat'] <-  dfm[k:(k+prod(Dim[-2])-1),'discards']/dfm[k:(k+prod(Dim[-2])-1),'catch']
            revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,6), sum, na.rm=T)[,years,]
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
mtStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + metier + stock + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
      res <- cbind(res[,1:6], data.frame(res[,7]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[7:(7+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + metier + stock + indicator + year + scenario + season, obj, quantile, prob = prob, na.rm=T)
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
                       quantile, prob = prob, na.rm=T)
      
      res <- cbind(res[,1:4], 
                   data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),  
                   data.frame(res[,8]),  data.frame(res[,9]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
      names(res)[-c(1:4)] <- unlist(mget(paste('nms', 1:5, sep="")))
    }
    else{
      res <- aggregate(list(catch = obj$catch, discards = obj$discards, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price),                   
                       list(fleet = obj$fleet, metier = obj$metier, stock = obj$stock, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.rm = TRUE)
      
      res <- cbind(res[,1:5], 
                   data.frame(res[,6]),  data.frame(res[,7]),   data.frame(res[,8]),  data.frame(res[,9]),
                   data.frame(res[,10]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
     
      names(res)[-c(1:5)] <-  unlist(mget(paste('nms', 1:5, sep="")))
      
    }
  }
  
  return(res)
}


#------------------------------------------------------------------------------#
# mtSum data.frame[year, season, stock, fleet, metier, iter ||,|| 
#        landings, discards, price] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
mtSum <- function(obj, flnms = names(obj$fleets),
                     years = dimnames(obj$biols[[1]]@n)[[2]], byyear = TRUE, long = TRUE, scenario = 'bc'){
  
  if(flnms[1] == 'all') flnms <- names(obj$fleets)
  
  fleets <- obj$fleets
  
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
                       income = numeric(n), 
                       vcost = numeric(n),
                       stringsAsFactors = FALSE)
    k <- 1
    for(m in mts){
      mt <- fl@metiers[[m]]
      dff[k:(k+prod(Dim)-1),'effort']   <- c((fl@effort*mt@effshare)[,years,])
      dff[k:(k+prod(Dim)-1),'effshare'] <- c(mt@effshare[,years,])
      dff[k:(k+prod(Dim)-1),'vcost']    <- c((fl@effort*mt@effshare*mt@vcost)[,years,])
       dff[k:(k+prod(Dim)-1),'income'] <- c(Reduce('+', lapply(mt@catches, function(x) unitSums(quantSums(x@landings.n*x@price))[,years])))

      k <- k + prod(Dim)
    }
    
    res <- rbind(res, dff)      
  }}
  else{
    for(f in flnms){
      fl <- fleets[[f]]
      mts <- names(fl@metiers)
      n <- prod(Dim)*length(mts)
      
      dff <-  data.frame(year = rep(years, prod(Dim[3])*length(mts)), 
                         fleet = rep(f, n), 
                         metier = rep(mts, each = prod(Dim)),
                         iter = rep(rep(1:Dim[3], each = prod(Dim[1])), length(mts)),  
                         effshare = numeric(n), 
                         effort = numeric(n),
                         income = numeric(n), 
                         vcost = numeric(n),
                         stringsAsFactors = FALSE)
      k <- 1
      for(m in mts){
        mt <- fl@metiers[[m]]
        dff[k:(k+prod(Dim[-2])-1),'effort']   <- c(seasonSums((fl@effort*mt@effshare)[,years,]))
        dff[k:(k+prod(Dim[-2])-1),'effshare'] <- c(seasonSums(mt@effshare[,years,]))
        dff[k:(k+prod(Dim[-2])-1),'vcost']    <- c(seasonSums(fl@effort*mt@effshare*mt@vcost)[,years,])
        dff[k:(k+prod(Dim[-2])-1),'income'] <- c(Reduce('+', lapply(mt@catches, function(x) seasonSums(unitSums(quantSums(x@landings.n*x@price)))[,years])))
        
        k <- k + prod(Dim)
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
mtSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 9){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + metier + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[6:(6+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + metier + indicator + year + scenario + season, obj, quantile, prob = prob, na.rm=T)
      res <- cbind(res[,1:6], data.frame(res[,7]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      names(res)[7:(7+length(prob)-1)] <- nms
    }
  }
  else{
    
    if(!('season' %in% names(obj))){
      res <- aggregate(list(effort = obj$effort,  effshare = obj$effshare, vcost = obj$vcost, income = obj$income,       
                            vcost = obj$vcost), 
                       list(fleet = obj$fleet, metier = obj$metier, year = obj$year, scenario = obj$scenario), 
                       quantile, prob = prob, na.rm=T)
      
      res <- cbind(res[,1:3], 
                   data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  
                   data.frame(res[,7]))
      
      nms1  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('effshare_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
       
      names(res)[-c(1:3)] <- unlist(mget(paste('nms', 1:4, sep="")))
    }
    else{
      res <- aggregate(list(effort = obj$effort,  effshare = obj$effshare, vcost = obj$vcost, income = obj$income,       
                            vcost = obj$vcost),                    
                       list(fleet = obj$fleet, metier = obj$metier,  year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.rm = TRUE)
      
      res <- cbind(res[,1:4], 
                   data.frame(res[,5]),  data.frame(res[,6]),   data.frame(res[,7]),  data.frame(res[,8]))
                   
      
      nms1  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('effshare_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:4)] <-  unlist(mget(paste('nms', 1:4, sep="")))
      
    }
  }
  
  return(res)
}

#------------------------------------------------------------------------------#
# advSum data.frame[year, season, fleet, metier, iter ||,|| 
#        effort, effshare] 
#------------------------------------------------------------------------------#
#' @rdname bioSum
advSum <- function(obj, stknms = catchNames(obj$fleets), 
                     years = dimnames(obj$biols[[1]]@n)[[2]], long = TRUE, scenario = 'bc'){
  
  
  fleets <- obj$fleets
  advice <- obj$advice
  
  if(stknms[1] == 'all') stknms <- names(biols)
  
  sts <- stknms
  
  Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,6)]
  Dimnm <- dimnames(fleets[[1]]@effort[,years,])
  
  res <- NULL
  
    
  n <- prod(Dim)*length(stknms)
        
  dfm <-  data.frame(year = rep(years, prod(Dim[2])*length(sts)), 
                           stock = rep(sts, each = prod(Dim)),
                           iter = rep(rep(1:Dim[2], each = prod(Dim[1])), length(sts)),  
                           catch = numeric(n), 
                           discards  = numeric(n),
                           discRat  = numeric(n),
                           landings = numeric(n), 
                           quotaUpt = numeric(n), 
                           tac = numeric(n),
                     stringsAsFactors = FALSE)
   k <- 1
        
   for(ss in sts){

      dfm[k:(k+prod(Dim)-1),'landings'] <- c(quantSums(seasonSums(unitSums(landWStock(fleets,ss)[,years]))))
      dfm[k:(k+prod(Dim)-1),'discards'] <- c(quantSums(seasonSums(unitSums(discWStock(fleets,ss)[,years]))))
      dfm[k:(k+prod(Dim)-1),'catch']    <- c(quantSums(seasonSums(unitSums(catchWStock(fleets,ss)[,years]))))
      dfm[k:(k+prod(Dim)-1),'discRat']  <-  dfm[k:(k+prod(Dim)-1),'discards']/dfm[k:(k+prod(Dim)-1),'catch']
      dfm[k:(k+prod(Dim)-1),'tac']      <- c(obj$advice$TAC[ss,years])
      dfm[k:(k+prod(Dim)-1),'quotaUpt'] <- dfm[k:(k+prod(Dim)-1),'catch']/dfm[k:(k+prod(Dim)-1),'tac']   
          
    k <- k + prod(Dim)
  }
  res <- rbind(res, dfm) 
       
  if(long == TRUE){ # transform res into long format
    r1 <-  4
    r2 <- 9
    
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
advSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 8){ # the object is in long format
    
    res <- aggregate(value ~  stock + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
      
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
    names(res)[5:(5+length(prob)-1)] <- nms
  }
  
  else{
    
      res <- aggregate(list(catch = obj$catch,  discards = obj$discards, discRat = obj$discRat, landings = obj$landings,       
                            quotaUpt = obj$quotaUpt, tac = obj$tac), 
                       list(stock = obj$stock, year = obj$year, scenario = obj$scenario), 
                       quantile, prob = prob, na.rm=T)
      
      res <- cbind(res[,1:3], 
                   data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  
                   data.frame(res[,7]),  data.frame(res[,8]),  data.frame(res[,9]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('tac_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:3)] <- unlist(mget(paste('nms', 1:6, sep="")))
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
riskSum <- function(obj, stknms = names(obj$biols), Bpa, Blim, Prflim, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], scenario = 'bc'){

  if (stknms == 'all') { 
    stknms <- names(obj$biols)
  } else if (sum(!stknms %in% names(obj$biols))>0) {
    stop(paste("'stknms' values should be in the following list: ", paste(names(obj$biols), collapse = ", "), sep=''))
  }
  if (sum(!names(Bpa) %in% stknms) + sum(!names(Blim) %in% stknms)>0) {
    stop(paste("Check names for 'Bpa' and 'Blim'. Values should be in the following list: ", paste(stknms, collapse = ", "), sep=''))
  }
  
  bioS <- bioSum(obj, stknms = stknms, years = dimnames(obj$biols[[1]]@n)[[2]], long = FALSE)
  bioS <- cbind(bioS, Bpa = Bpa[bioS$stock],  Blim = Blim[bioS$stock])
  bioS <- cbind(bioS, risk.pa = as.numeric(bioS$ssb < bioS$Bpa), risk.lim = as.numeric(bioS$ssb < bioS$Blim))

  if (flnms == 'all') { 
    flnms <- names(obj$fleets)
  } else if (sum(!flnms %in% names(obj$fleets))>0) {
    stop(paste("'flnms' values should be in the following list: ", paste(names(obj$fleets), collapse = ", "), sep=''))
  }
  if (sum(!names(Prflim) %in% flnms)>0) {
    stop(paste("Check names for 'Prflim'. Values should be in the following list: ", paste(flnms, collapse = ", "), sep=''))
  }
  
  flS <- fltSum(obj, years = dimnames(obj$biols[[1]]@n)[[2]], flnms = flnms, long = FALSE)
  flS <- cbind(flS, refp = Prflim[flS$fleet])
  flS <- cbind(flS, risk = as.numeric(flS$profits < flS$refp))
  
  if(all(is.na(flS$risk))){ # if economic data not available for example
    flS$risk <- 0
    auxFl    <- aggregate(risk ~ year + fleet + scenario, data=flS, FUN=function(x){sum(x)/length(x)})
    auxFl$risk[] <- NA
  }else{
    auxFl    <- aggregate(risk ~ year + fleet + scenario, data=flS, FUN=function(x){sum(x)/length(x)})
  }
  
  # auxBioPa   <- aggregate(risk.pa ~ year + stock + scenario, data=bioS, FUN=function(x){sum(x)/length(x)})
  if(all(is.na(bioS$risk.pa))){ # if Bpa not available
    bioS$risk.pa <- 0
    auxBioPa     <- aggregate(risk.pa ~ year + stock + scenario, data=bioS, FUN=function(x){sum(x)/length(x)})
    auxBioPa$risk.pa[] <- NA
  }else{
    auxBioPa     <- aggregate(risk.pa ~ year + stock + scenario, data=bioS, FUN=function(x){sum(x)/length(x)})
  }
  
  # auxBiolim  <- aggregate(risk.lim ~ year + stock + scenario, data=bioS, FUN=function(x){sum(x)/length(x)})
  if(all(is.na(bioS$risk.lim))){ # if Blim not available
    bioS$risk.lim <- 0
    auxBiolim     <- aggregate(risk.lim ~ year + stock + scenario, data=bioS, FUN=function(x){sum(x)/length(x)})
    auxBiolim$risk.lim[] <- NA
  }else{
    auxBiolim     <- aggregate(risk.lim ~ year + stock + scenario, data=bioS, FUN=function(x){sum(x)/length(x)})
  }
  
  names(auxFl) <- c('year', 'unit', 'scenario', 'value')
  names(auxBioPa) <- c('year', 'unit', 'scenario', 'value')
  names(auxBiolim) <- c('year', 'unit', 'scenario', 'value')
  
  res <- cbind(scenario = scenario, rbind(
               cbind(auxFl[,1:3],     indicator = 'pPrflim', value = auxFl[,4]),
               cbind(auxBioPa[,1:3],  indicator = 'pBpa',    value = auxBioPa[,4]),
               cbind(auxBiolim[,1:3], indicator = 'pBlim',   value = auxBiolim[,4])))
 # No sense in wide format  
#  if(long == FALSE){
#    temp <- reshape(res,v.names = 'value', timevar = 'indicator', idvar = c('scenario', 'year', 'unit'), direction = 'wide')  }

  return(res)
}

#----------------------------------------------------------------------
# npv(obj, years, flnms)
#----------------------------------------------------------------------
#' @rdname bioSum
npv <- function(obj, discF = 0.05, y0, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], scenario = 'bc'){
  
  flS <- fltSum(obj, years = years, flnms = flnms, long = FALSE)
  
  y0 <- as.numeric(y0)
  
  flS <- cbind(flS, discount= (1+discF)^(as.numeric(as.character(flS$year))-y0))
  flS <- cbind(flS, discProf = flS$profits/(flS$discount))
  
  flS <- subset(flS, year %in% c(years))
  
  res <- aggregate(discProf ~ fleet + iter, data = flS, FUN = sum)
  
  names(res)[3] <- 'npv'
  
  res <- cbind(scenario = scenario, res)
  
  return(res)
  
}  
    

#' @rdname bioSum
npvQ <- function(obj, prob = c(0.05,0.5,0.95)){
  
  res <- aggregate(npv ~  fleet + scenario, obj, quantile, prob = prob,na.rm=T)
  res <- cbind(fleet = res[,1], data.frame(res[,2]))
  
  nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
  
  names(res)[2:(2+length(prob)-1)] <- nms
  
  
  return(res)
  
}   


#----------------------------------------------------------------------
# vesselSum
#----------------------------------------------------------------------
#' @rdname bioSum
vesselSum <- function(obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = TRUE, scenario = 'bc'){
  
  flS <- fltSum(obj, flnms = flnms, years = years, byyear = byyear, long = long, scenario = scenario)
  
  if(byyear == TRUE){
  if(long == FALSE){
    for(col in 5:21){ 
      flS[,col] <- flS[,col]/flS[,'nVessels']
    }
    res <- flS[,c(1:4,6:15,17:22)]
  }
  else{
   ids <- c("catch", "costs", "discards","discRat","effort",  "fcosts","gva","income", "landings",     
      "netProfit",  "price", "profits", "quotaUpt", "salaries", "vcosts")
   for(id in ids){
     flS[flS$indicator == id, 'value'] <- flS[flS$indicator == id, 'value']/flS[flS$indicator == 'nVessels', 'value']
   }
   res <- flS
  }}
   else{
     if(long == FALSE){
       for(col in 6:22){ 
         flS[,col] <- flS[,col]/flS[,'nVessels']
       }
       res <- flS[,c(1:5,7:16,18:23)]
     }
     else{
       ids <- c("catch", "costs", "discards","discRat","effort",  "fcosts","gva","income", "landings",     
                "netProfit",  "price", "profits", "quotaUpt", "salaries", "vcosts")
       for(id in ids){
         flS[flS$indicator == id, 'value'] <- flS[flS$indicator == id, 'value']/flS[flS$indicator == 'nVessels', 'value']
       }
       res <- subset(flS, indicator %in% c(ids, 'profitability'))
   }
   
  }
     return(res) 
               
}

#' @rdname bioSum
vesselSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){

if(dim(obj)[2] < 10){ # the object is in long format
  
  if(!('season' %in% names(obj))){
    res <- aggregate(value ~ fleet + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
    names(res)[5:(5+length(prob)-1)] <- nms
  }
  else{
    res <- aggregate(value ~ fleet + indicator + year + scenario + season, obj, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:5], data.frame(res[,6]))
    
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[6:(6+length(prob)-1)] <- nms
  }
}
else{
  
  if(!('season' %in% names(obj))){
    res <- aggregate(list(    catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                          discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                          income  = obj$income,         landings  = obj$landings,  netProfit  = obj$netProfit,    
                          price  = obj$price,           profits  = obj$profits,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                          vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                     list(fleet = obj$fleet, year = obj$year, scenario = obj$scenario), 
                     quantile, prob = prob, na.rm=T)
    
    res <- cbind(res[,1:3], 
                 data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                 data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                 data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                 data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]))
    
     nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms3  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms4  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms5  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms6  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms7  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms8  <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms9  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms10 <- paste('netprofit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms12 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms13 <- paste('profit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms14 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms15 <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms16 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms17 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms18 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
    names(res)[-c(1:3)] <- unlist(mget(paste('nms', c(2:10,12:18), sep="")))
  }
  else{
    res <- aggregate(list(    catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                          discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                          income  = obj$income,         landings  = obj$landings,  netProfit  = obj$netProfit,    
                          price  = obj$price,           profits  = obj$profits,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                          vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                     list(fleet = obj$fleet, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.rm = TRUE)
    
    res <- cbind(res[,1:3], 
                 data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                 data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                 data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                 data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]))
    
    nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms3  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms4  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms5  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms6  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms7  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms8  <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms9  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms10 <- paste('netprofit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms12 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms13 <- paste('profit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms14 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms15 <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms16 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms17 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    nms18 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    
    names(res)[-c(1:3)] <- unlist(mget(paste('nms',  c(2:10,12:18), sep="")))
    
  }
}

return(res)
}


#----------------------------------------------------------------------
# vesselSum
#----------------------------------------------------------------------
#' @rdname bioSum
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
      for(col in c(6:9,11,13)){
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
      for(col in c(6:9,11,13)){
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
vesselStkSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 10){ # the object is in long format
    
    if(!('season' %in% names(obj))){
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario, obj, quantile, prob = prob,na.rm=T)
      res <- cbind(res[,1:5], data.frame(res[,6]))
      
      nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[6:(6+length(prob)-1)] <- nms
    }
    else{
      res <- aggregate(value ~ fleet + stock + indicator + year + scenario + season, obj, quantile, prob = prob, na.rm=T)
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
                       quantile, prob = prob, na.rm=T)
      
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
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, season = obj$season, scenario = obj$scenario), quantile, prob = prob, na.rm = TRUE)
      
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


