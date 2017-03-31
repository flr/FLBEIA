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

#' Biological summary functions
#' 
#' These functions return the biomass (B), fishing mortality (F),  spawning stock biomass (SSB), recruitment (R), catches (C), landings (L) and discards (D) indicators. 
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
      # spawning season: first season with fraction of natural mortality before spawning < 1
      spwn.sson <- 1
      si <- 0
      while( (si-spwn.sson)!=0) { 
        si <- spwn.sson
        spwn.sson  <- ifelse( sum(spwn(obj$biols[[stk]])[ , , 1, spwn.sson, drop = T]<1,na.rm=T)==0, spwn.sson+1, spwn.sson)
        d  <- si-spwn.sson 
      }
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
#      \item{Summary_flbeia:} The data frame contains five columns: stock, years, iteration, 
#      indicator and value. The indicators are: recruitment, ssb, f, biomass, catch, landings and discards.
#'      \item{bioSum:} Data frame with the biological indicators The columns are stocks, 
#'      years, it, indicators and value. 
#'      \item{fleetSum:} The colums of the data frame are: year, season, stock, fleet, iter, profits, 
#'      capacity, costs, discards, effort and landings.
#'      \item{ecoSum_damara:} ecoSum for Damara project.
#'      \item{metierSum:} The colums of the data frame are: year, season, fleet, metier, iter, 
#'      effort and effshare.
#'      \item{stockFleetSum:} The colums of the data frame are: year, season, stock, fleet, iter, 
#'      landings, discards, price and tacshare.
#'      \item{stockMetierSum:} The colums of the data frame are: year, season, stock, fleet, metier, iter, 
#'      landings, discards and price.
#'      
#'}
#'
#' @return The data frames can be of wide or long format. In long format all the indicators are in the same column.
#'  There is one column, indicator, for the name of the indicator and a second one value for the numerica value of the indicator.
#'  In the wide format each of the indicators correspond with one column in the data frame. 
#'  The long format it is recommendable to work with ggplot2 functions for example and the wide format 
#'  it is more efficient for memory allocation and speed of computations. 

#' @arguments
#' @inheritParams FLBEIA
#' @inheritParams F_flbeia
#' @param flnms Names of the fleets.
#' @param stknms Names of the stocks.
#' @param years The period of years. 
#' @param long The data frame should be constructed using long or wide format?
#' 
#
#' @return A data frame.

#' @examples
#'\dontrun{
#'
#' library(FLBEIA)
#' library(ggplot2)
#'
#' # Apply the summary functions to the examples runs in FLBEIA help page.
#' # Test the different arguments in summary function.
#' 
#' data(res_flbeia)
#'
#' # Example One: One stock, one fleet, one iter.
#' s0_bio    <- bioSum(s0)
#' s0_flt    <- fltSum(s0)
#' 
#' s0_bioQ   <- bioSumQ(s0_bio)
#' s0_fltQ   <- fltSumQ(s0_flt)
#' 
#' s0_bio    <- bioSum(s0, long = FALSE, years = ac(2016:2020))
#' s0_flt    <- fltSum(s0, long = FALSE, years = ac(2016:2020))
#' 
#' s0_bioQ   <- bioSumQ(s0_bio)
#' s0_fltQ   <- fltSumQ(s0_flt)
#' 
#' s0_flt    <- fltSum(s0, long = FALSE, byyear = FALSE)
#' s0_fltQ   <- fltSumQ(s0_flt)
#' 
#' 
#' # Example OneIters: As one but with iterations.
#' s1_bio    <- bioSum(s1, years = ac(2016:2020))
#' s1_flt    <- fltSum(s1, years = ac(2016:2020))
#' 
#' s1_bioQ   <- bioSumQ(s1_bio)
#' s1_fltQ   <- fltSumQ(s1_flt)
#' 
#' s1_bio    <- bioSum(s1, long = FALSE)
#' s1_flt    <- fltSum(s1, long = FALSE)
#' 
#' s1_bioQ   <- bioSumQ(s1_bio)
#' s1_fltQ   <- fltSumQ(s1_flt)
#' 
#' s1_flt    <- fltSum(s1, long = FALSE, byyear = FALSE)
#' s1_fltQ   <- fltSumQ(s1_flt)
#' 
#' # Example Multi: Two stock, two fleet, four iters. 
#' s2_bio    <- bioSum(s2, years = ac(2016:2020))
#' s2_flt    <- fltSum(s2, years = ac(2016:2020))
#' 
#' s2_bioQ   <- bioSumQ(s2_bio)
#' s2_fltQ   <- fltSumQ(s2_flt)
#' 
#' s2_bio    <- bioSum(s2, long = FALSE)
#' s2_flt    <- fltSum(s2, long = FALSE)
#' 
#' s2_bioQ   <- bioSumQ(s2_bio)
#' s2_fltQ   <- fltSumQ(s2_flt)
#' 
#' s2_flt    <- fltSum(s2, long = FALSE, byyear = FALSE)
#' s2_fltQ   <- fltSumQ(s2_flt)
#' 
#' s2_bio    <- bioSum(s2, stocks = 'stk2')
#' s2_flt    <- fltSum(s2, flnms = 'fl1')
#' 
#' s2_bioQ   <- bioSumQ(s2_bio)
#' s2_fltQ   <- fltSumQ(s2_flt)
#' 
#' }
bioSum <- function(obj, stocks = 'all', years = dimnames(obj$biols[[1]]@n)$year, long = TRUE, scenario = 'bc'){
    xx <- summary_flbeia(obj, years)
     
    dnms <- dimnames(xx)
    
    if(stocks == 'all') stocks <- dnms[[1]]
    
    xx <- xx[stocks,,,,drop=F]
    dnms <- dimnames(xx)
    
    if(long == TRUE){
      
      df <- expand.grid(iter = dnms[[3]], year = dnms[[2]],indicator = c(dnms[[4]], 'land.iyv', 'disc.iyv', 'catch.iyv'),  stock = stocks)[,4:1]

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
        }
      }
    }
    
    df <- cbind(scenario = scenario, df)
    return(df)
}

bioSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){

  if(dim(obj)[2] <= 6){ # the object is in long format
    res <- aggregate(value ~ stock + indicator + year + scenario, obj, quantile, prob = prob)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
           
    names(res)[5:(5+length(prob)-1)] <- nms
  }
  else{
    res <- aggregate(list(biomass = obj$biomass,catch = obj$catch,catch.iyv = obj$catch.iyv, 
                          discards = obj$discards, disc.iyv = obj$disc.iyv,
                          f = obj$f,landings = obj$landings, land.iyv = obj$land.iyv,
                          rec = obj$rec, ssb = obj$ssb), 
                       list(stock = obj$stock, year = obj$year), quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]), data.frame(res[,4]), data.frame(res[,5]),
                           data.frame(res[,6]), data.frame(res[,7]), data.frame(res[,8]),
                           data.frame(res[,9]), data.frame(res[,10]), data.frame(res[,11]),
                           data.frame(res[,12]))
                                                                                 
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
    
    names(res)[-c(1:2)] <- c(nms_bio, nms_cat, nms_iyvcat, nms_disc, nms_iyvdisc, nms_f, nms_land, nms_iyvland,nms_rec, nms_ssb)
    
    
  }
  
  return(res)
  }
  
  
#------------------------------------------------------------------------------#
# ecoSum data.frame[year, season, stock, fleet, iter, ||,|| 
#        profits, capacity, costs, discards, effort, landings] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia
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
      covars$MaxDays <- FLQuant(0, dim = c(length(fleets),dim(fleets[[1]]@effort[, years, ])[2:6]))
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
        vcosts   = numeric(n))
    
    k <- 1
    for (f in flnms) {
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(catchWStock.f(fl, x))))
        res[k:(k + prod(Dim) - 1), "catch"] <- c(Reduce('+',temp))
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(landWStock.f(fl, x))))
        res[k:(k + prod(Dim) - 1), "landings"] <- c(Reduce('+',temp))
        
        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(discWStock.f(fl, x))))
        res[k:(k + prod(Dim) - 1), "discards"] <- c(Reduce('+',temp))
        
        res[k:(k + prod(Dim) - 1), "discRat"] <- res[k:(k + prod(Dim) - 1), "discards"]/res[k:(k + prod(Dim) - 1), "catch"]
          
        res[k:(k + prod(Dim) - 1), "capacity"] <- c(fl@capacity[,years, ])
        
        res[k:(k + prod(Dim) - 1), "effort"] <- c(fl@effort[,years, ])
        
        res[k:(k + prod(Dim) - 1), "fcosts"] <- c(totfcost_flbeia(fl, covars, f)[,years, ])
        
        res[k:(k + prod(Dim) - 1), "vcosts"] <- c(totvcost_flbeia(fl)[,years, ])
        
        res[k:(k + prod(Dim) - 1), "costs"] <- c(costs_flbeia(fl, covars, f)[,years, ])
       
        res[k:(k + prod(Dim) - 1), "income"] <- c(revenue_flbeia(fl)[,years, ]) 
        
        res[k:(k + prod(Dim) - 1), "profits"] <- c(revenue_flbeia(fl)[,years, ]) - res[k:(k + prod(Dim) - 1), "costs"]
         
        res[k:(k + prod(Dim) - 1), "salaries"] <- c(fl@crewshare[,years,]*revenue_flbeia(fl)[,years, ] + covars[['Salaries']][f,years])
        
        res[k:(k + prod(Dim) - 1), "gva"] <- res[k:(k + prod(Dim) - 1), "income"] - res[k:(k + prod(Dim) - 1), "costs"] + res[k:(k + prod(Dim) - 1), "salaries"]
        
        res[k:(k + prod(Dim) - 1), "profitability"] <- res[k:(k + prod(Dim) - 1), "profits"]/res[k:(k + prod(Dim) - 1), "income"]
 
        res[k:(k + prod(Dim) - 1), "nVessels"]  <- c(covars[['NumbVessels']][f,])
          
        res[k:(k + prod(Dim) - 1), "netProfit"] <- c(revenue_flbeia(fl)[,years, ] -  costs_flbeia(fl, covars, f)[,years, ] - covars[['Depreciation']][f,]*covars[['NumbVessels']][f,])

        temp <- lapply(catchNames(fl), function(x) quantSums(unitSums(catchWStock.f(fl, x))))
        temp <- Reduce('+',temp)
        totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,]*obj$advice$TAC[x,]))
        res[k:(k + prod(Dim) - 1), "quotaUpt"] <- c(sweep(temp, 4,totTAC/dim(temp)[4], "/"))
        
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
                        vcosts   = numeric(n))
      
      k <- 1
      for (f in flnms) {
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
        res[k:(k + prod(Dim) - 1), "catch"] <- c(Reduce('+',temp)[, years, ])
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(landWStock.f(fl, x)))))
        res[k:(k + prod(Dim) - 1), "landings"] <- c(Reduce('+',temp)[, years, ])
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(discWStock.f(fl, x)))))
        res[k:(k + prod(Dim) - 1), "discards"] <- c(Reduce('+',temp)[, years, ])
        
        res[k:(k + prod(Dim) - 1), "discRat"] <- res[k:(k + prod(Dim) - 1), "discards"]/res[k:(k + prod(Dim) - 1), "catch"]
        
        res[k:(k + prod(Dim) - 1), "capacity"] <- c(seasonSums(fl@capacity[,years, ]))
        
        res[k:(k + prod(Dim) - 1), "effort"] <- c(seasonSums(fl@effort[,years, ]))
        
        res[k:(k + prod(Dim) - 1), "fcosts"] <- c(seasonSums(totfcost_flbeia(fl, covars, f)[,years, ]))
        
        res[k:(k + prod(Dim) - 1), "vcosts"] <- c(seasonSums(totvcost_flbeia(fl)[,years, ]))
        
        res[k:(k + prod(Dim) - 1), "costs"] <- c(seasonSums(costs_flbeia(fl, covars, f)[,years, ]))
        
        res[k:(k + prod(Dim) - 1), "income"] <- c(seasonSums(revenue_flbeia(fl)[,years, ])) 
        
        res[k:(k + prod(Dim) - 1), "profits"] <- c(seasonSums(revenue_flbeia(fl)[,years, ])) - res[k:(k + prod(Dim) - 1), "costs"]
        
        res[k:(k + prod(Dim) - 1), "salaries"] <- c(seasonSums(fl@crewshare[,years,]*revenue_flbeia(fl)[,years, ] + covars[['Salaries']][f,years]))
        
        res[k:(k + prod(Dim) - 1), "gva"] <- res[k:(k + prod(Dim) - 1), "income"] - res[k:(k + prod(Dim) - 1), "costs"] + res[k:(k + prod(Dim) - 1), "salaries"]
        
        res[k:(k + prod(Dim) - 1), "profitability"] <- res[k:(k + prod(Dim) - 1), "profits"]/res[k:(k + prod(Dim) - 1), "income"]
        
        res[k:(k + prod(Dim) - 1), "nVessels"]  <- c(seasonMeans(covars[['NumbVessels']][f, years, ]))
        
        res[k:(k + prod(Dim) - 1), "netProfit"] <- c(seasonSums(revenue_flbeia(fl)[,years, ] -  costs_flbeia(fl, covars, f)[,years, ] - covars[['Depreciation']][f,years]*covars[['NumbVessels']][f,years]))
        
        temp <- lapply(catchNames(fl), function(x) seasonSums(quantSums(unitSums(catchWStock.f(fl, x)))))
        temp <- Reduce('+',temp)[, years, ]
        totTAC <- Reduce('+',lapply(names(obj$advice$quota.share), function(x) obj$advice$quota.share[[x]][f,years]*obj$advice$TAC[x,years]))
        res[k:(k + prod(Dim) - 1), "quotaUpt"] <- c(totTAC/temp)
        
        k <- k + prod(Dim)
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
                          list(fleet = obj$fleet, year = obj$year), 
                          quantile, prob = prob, na.rm=T)
    
      res <- cbind(res[,1:2], 
                 data.frame(res[,3]),  data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),
                 data.frame(res[,7]),  data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),
                 data.frame(res[,11]), data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]),
                 data.frame(res[,15]), data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]),
                 data.frame(res[,19]), data.frame(res[,20]))
                 
      nms1  <- paste('capacity_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms8  <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms9  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms10 <- paste('netprofit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms11 <- paste('nVessels_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms12 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms13 <- paste('profit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms14 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms15 <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms16 <- paste('salaries_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms17 <- paste('vcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms18 <- paste('profitability_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
   
      names(res)[-c(1:2)] <- unlist(mget(paste('nms', 1:18, sep="")))
    }
    else{
        res <- aggregate(list(capacity = obj$capacity,      catch = obj$catch,         costs = obj$costs,          discards = obj$discards,       
                              discRat = obj$discRat,    effort = obj$effort,       fcosts = obj$fcosts,        gva  = obj$gva,                    
                              income  = obj$income,         landings  = obj$landings,  netProfit  = obj$netProfit, nVessels  = obj$nVessels,     
                              price  = obj$price,           profits  = obj$profits,    quotaUpt  = obj$quotaUpt,   salaries  = obj$salaries, 
                              vcosts  = obj$vcosts,         profitability  = obj$profitability), 
                         list(fleet = obj$fleet, year = obj$year, season = obj$season), quantile, prob = prob, na.rm = TRUE)
        
        res <- cbind(res[,1:3], 
                     data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                     data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]),  data.frame(res[,11]),
                     data.frame(res[,12]), data.frame(res[,13]), data.frame(res[,14]), data.frame(res[,15]),
                     data.frame(res[,16]), data.frame(res[,17]), data.frame(res[,18]), data.frame(res[,19]),
                     data.frame(res[,20]), data.frame(res[,21]))
        
        nms1  <- paste('capacity_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms2  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms3  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms4  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms5  <- paste('effort_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms6  <- paste('fcost_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms7  <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms8  <- paste('income_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms9  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms10 <- paste('netprofit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms11 <- paste('nVessels_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms12 <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms13 <- paste('profit_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms14 <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
        nms15 <- paste('gva_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
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
     return(sweep(fleet@fcost, 4,covars[["NumbVessels"]][flnm, ], "*"))            
}


#------------------------------------------------------------------------------#
# catchFlSum data.frame[year, season, stock, fleet, iter, ||,|| 
#        landings, discards, price, tacshare] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia
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
                    tacshare = numeric(n))
        
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
                          tacshare = numeric(n))
        
        k <- 1
        
        for(st in sts){
          
          dff[k:(prod(Dim) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,6), sum)[,years])    
          dff[k:(prod(Dim) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,6), sum)[,years]) 
          dff[k:(prod(Dim) + k-1),'catch']    <- dff[k:(prod(Dim) + k-1),'discards'] + dff[k:(prod(Dim) + k-1),'landings']
          dff[k:(prod(Dim) + k-1),'discRat']  <- dff[k:(prod(Dim) + k-1),'discards']/dff[k:(prod(Dim) + k-1),'catch']
          dff[k:(prod(Dim) + k-1),'price']    <- c(seasonMeans(price_flbeia(fl, st)[,years]*quantSums(landWStock.f(fl, st)[,years]))/seasonSums(quantSums(landWStock.f(fl, st)[,years])))
          dff[k:(prod(Dim) + k-1),'quota']    <- c((advice$TAC[st,]*advice$quota.share[[st]][f,])[,years])
          dff[k:(prod(Dim) + k-1),'quotaUpt'] <- dff[k:(prod(Dim) + k-1),'catch']/dff[k:(prod(Dim) + k-1),'quota']
          
          k <- k + prod(Dim)     
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
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year), 
                       quantile, prob = prob, na.rm=T)
      
      res <- cbind(res[,1:3], 
                   data.frame(res[,4]),  data.frame(res[,5]),  data.frame(res[,6]),  data.frame(res[,7]),
                   data.frame(res[,8]),  data.frame(res[,9]),  data.frame(res[,10]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('price_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('quota_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms7  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
     
      names(res)[-c(1:3)] <- unlist(mget(paste('nms', 1:7, sep="")))
    }
    else{
      res <- aggregate(list(catch = obj$catch, discards = obj$catch, discRat = obj$discRat, landings = obj$landings,       
                            price = obj$price, quota = obj$quota,    quotaUpt = obj$quotaUpt),                   
                       list(fleet = obj$fleet, stock = obj$stock, year = obj$year, season = obj$season), quantile, prob = prob, na.rm = TRUE)
      
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
#' @rdname summary_flbeia 
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
                        price = numeric(n))
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
                             price = numeric(n))
          k <- 1
          
          for(ss in sts){
            cc <- mt@catches[[ss]]
            dfm[k:(k+prod(Dim)-1),'landings'] <- c(apply(cc@landings[,years,], c(2,6), sum, na.rm=T))
            dfm[k:(k+prod(Dim)-1),'discards'] <- c(apply(cc@discards[,years,], c(2,6), sum, na.rm=T))
            dfm[k:(k+prod(Dim)-1),'catch'] <- dfm[k:(k+prod(Dim)-1),'discards'] + dfm[k:(k+prod(Dim)-1),'landings']
            dfm[k:(k+prod(Dim)-1),'discRat'] <-  dfm[k:(k+prod(Dim)-1),'discards']/dfm[k:(k+prod(Dim)-1),'catch']
            revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,6), sum, na.rm=T)[,years,]
            dfm[k:(k+prod(Dim)-1),'price']  <- c(revst)/dfm[k:(k+prod(Dim)-1),'landings']  
            k <- k + prod(Dim)
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
                       list(fleet = obj$fleet, metier = obj$metier, stock = obj$stock, year = obj$year), 
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
                       list(fleet = obj$fleet, metier = obj$metier, stock = obj$stock, year = obj$year, season = obj$season), quantile, prob = prob, na.rm = TRUE)
      
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
#' @rdname summary_flbeia 
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
                       vcost = numeric(n))
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
                         vcost = numeric(n))
      k <- 1
      for(m in mts){
        mt <- fl@metiers[[m]]
        dff[k:(k+prod(Dim)-1),'effort']   <- c(seasonSums((fl@effort*mt@effshare)[,years,]))
        dff[k:(k+prod(Dim)-1),'effshare'] <- c(seasonSums(mt@effshare[,years,]))
        dff[k:(k+prod(Dim)-1),'vcost']    <- c(seasonSums(fl@effort*mt@effshare*mt@vcost)[,years,])
        dff[k:(k+prod(Dim)-1),'income'] <- c(Reduce('+', lapply(mt@catches, function(x) seasonSums(unitSums(quantSums(x@landings.n*x@price)))[,years])))
        
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
mtSumQ <- function(obj,  prob = c(0.95,0.5,0.05)){
  
  if(dim(obj)[2] < 8){ # the object is in long format
    
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
                       list(fleet = obj$fleet, metier = obj$metier, year = obj$year), 
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
                       list(fleet = obj$fleet, metier = obj$metier,  year = obj$year, season = obj$season), quantile, prob = prob, na.rm = TRUE)
      
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
#' @rdname summary_flbeia
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
                           tac = numeric(n))
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
                       list(stock = obj$stock, year = obj$year), 
                       quantile, prob = prob, na.rm=T)
      
      res <- cbind(res[,1:2], 
                   data.frame(res[,3]),  data.frame(res[,4]),  data.frame(res[,5]),  
                   data.frame(res[,6]),  data.frame(res[,7]),  data.frame(res[,8]))
      
      nms1  <- paste('catch_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms2  <- paste('discards_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms3  <- paste('discRat_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms4  <- paste('landings_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms5  <- paste('quotaUpt_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      nms6  <- paste('tac_q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
      
      names(res)[-c(1:2)] <- unlist(mget(paste('nms', 1:6, sep="")))
    }
  
  return(res)
}

#----------------------------------------------------------------------
# riskSum(obj, stocks, fleets, years, long)
# Bpa = a named vector with the precautionary biomass per stock.
# Blim = a named vector with the limit biomass per stock.
# Blim = a named vector with the limit profit per fleet
#----------------------------------------------------------------------

riskSum <- function(obj, stocks = names(obj$biols), Bpa, Blim, Prflim, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], long = TRUE, scenario = 'bc'){

  bioS <- bioSum(obj, stocks = names(obj$biols), years = dimnames(obj$biols[[1]]@n)[[2]], long = FALSE)
  bioS <- cbind(bioS, Bpa = Bpa[bioS$stock],  Blim = Blim[bioS$stock])
  bioS <- cbind(bioS, risk.pa = as.numeric(bioS$ssb > bioS$Bpa), risk.lim = as.numeric(bioS$ssb > bioS$Blim))
  
  flS <- fltSum(obj, years = dimnames(obj$biols[[1]]@n)[[2]], flnms = names(obj$fleets), long = FALSE)
  flS <- cbind(flS, refp = Prflim[flS$fleet])
  flS <- cbind(flS, risk = as.numeric(flS$profits > flS$refp))
  
  auxFl      <- aggregate(risk ~  year + fleet, data=flS, FUN=function(x){sum(x)/length(x)})
  auxBioPa   <- aggregate(risk.pa ~ year + stock, data=bioS, FUN=function(x){sum(x)/length(x)})
  auxBiolim  <- aggregate(risk.lim ~ year + stock, data=bioS, FUN=function(x){sum(x)/length(x)})
  
  names(auxFl) <- c('year', 'unit', 'value')
  names(auxBioPa) <- c('year', 'unit', 'value')
  names(auxBiolim) <- c('year', 'unit', 'value')
  
  res <- cbind(scenario = scenario, rbind(
               cbind(auxFl[,1:2],     indicator = 'pPrflim', value = auxFl[,3]),
               cbind(auxBioPa[,1:2],  indicator = 'pBpa',    value = auxBioPa[,3]),
               cbind(auxBiolim[,1:2], indicator = 'pBlim',   value = auxBiolim[,3])))
  return(res)
}

#----------------------------------------------------------------------
# npv(obj, years, flnms)
#----------------------------------------------------------------------
npv <- function(obj, discF = 0.05, y0, flnms = names(obj$fleets), years = dimnames(obj$biols[[1]]@n)[[2]], scenario = 'bc'){
  
  flS <- fltSum(obj, years = dimnames(obj$biols[[1]]@n)[[2]], flnms = names(obj$fleets), long = FALSE)
  
  flS <- cbind(flS, discount= (1+discF)^(as.numeric(flS$year)-y0))
  flS <- cbind(flS, discProf = flS$profits/flS$discount)
  
  flS <- subset(flS, year %in% c(years))
  
  res <- aggregate(discProf ~ fleet, data = flS, FUN = sum)
  
  names(res)[2] <- 'npv'
  
  return(res)
  
}  
           