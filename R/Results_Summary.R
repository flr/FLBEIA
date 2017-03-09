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
#' These functions return the summary of the biomass (B), fishing mortality (F),  spawning stock biomass (SSB), recruitment (R), catches (C), landings (L) and discards (D). 
#' 
#' @param  obj The output of the FLBEIA function.
#' @param  years Is the period of years defined in the obj.
#' 
#' @return The summary of the selected biological item.
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
#'      \item{Summary_flbeia:} The data frame contains five columns: stock, years, iteration, 
#'      indicator and value. The indicators are: recruitment, ssb, f, biomass, catch, landings and discards.
#'      \item{bioSum:} The data frame contains biological results. The columns are stocks, 
#'      years, it, indicators and value. 
#'      \item{ecoSum:} The colums of the data frame are: year, quarter, stock, fleet, iter, profits, 
#'      capacity, costs, discards, effort and landings.
#'      \item{ecoSum_damara:} ecoSum for Damara project.
#'      \item{effortMtSum:} The colums of the data frame are: year, quarter, fleet, metier, iter, 
#'      effort and effshare.
#'      \item{catchFlSum:} The colums of the data frame are: year, quarter, stock, fleet, iter, 
#'      landings, discards, price and tacshare.
#'      \item{catchMtSum:} The colums of the data frame are: year, quarter, stock, fleet, metier, iter, 
#'      landings, discards and price.
#'      
#'}

#' @arguments
#' @inheritParams FLBEIA
#' @inheritParams F_flbeia
#' @param flnm Names of the fleets.
#' @param stknms Names of the stocks.
#' @param years The period of years. 
#' 
#
#' @return A data frame.

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(ggplot2)
#' data(one)
#' s0 <- FLBEIA(biols = oneBio,       # FLBiols object with one FLBiol element for stk1.
#'                SRs = oneSR,        # A list with one FLSRSim object for stk1.
#'                BDs = NULL,         # No Biomass Dynamic populations in this case.
#'             fleets = oneFl,        # FLFleets object with on fleet.
#'             covars = NULL,         # covars not used
#'            indices = NULL,         # indices not used 
#'             advice = oneAdv,       # A list with two elements 'TAC' and 'quota.share'
#'          main.ctrl = oneMainC,     # A list with one element to define the start and end of the simulation.
#'         biols.ctrl = oneBioC,      # A list with one element to select the model to simulate the stock dynamics.
#'        fleets.ctrl = oneFlC,       # A list with several elements to select fleet dynamic models and store additional parameters.
#'        covars.ctrl = NULL,         # covars control not used 
#'           obs.ctrl = oneObsC,      # A list with one element to define how the stock observed ("PerfectObs").
#'        assess.ctrl = oneAssC,      # A list with one element to define how the stock assessment model used ("NoAssessment").
#'        advice.ctrl = oneAdvC) 
#' summary_flbeia(obj=s0)
#' }
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
#' @rdname summary_flbeia
bioSum <- function(obj, years = dimnames(obj$biols[[1]]@n)$year){
    xx <- summary_flbeia(obj, years)
    n  <- prod(dim(xx))
    
    dnms <- dimnames(xx)
    
    df <- expand.grid(iter = dnms[[3]], indicator = dnms[[4]], year = dnms[[2]], stock = dnms[[1]])[,4:1]

    df$stock     <- as.character(df$stock)
    df$year      <- as.numeric(as.character(df$year))
    df$indicator <- as.character(df$indicator)
    df$iter      <- as.numeric(df$iter)

    df <- cbind(df, value = NA)
    
    for(st in dnms[[1]]){
        for(yr in dnms[[2]]){
            for(ind in dnms[[4]]){
               df[df$stock == st & df$year == yr & df$indicator == ind,'value'] <- xx[st,yr,,ind]

            }
        }
    }
    
    return(df)
}

#------------------------------------------------------------------------------#
# ecoSum data.frame[year, quarter, stock, fleet, iter, ||,|| 
#        profits, capacity, costs, discards, effort, landings] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia
ecoSum <- function (fleets, flnms = "all", years, covars = NULL)
{
    if (flnms[1] == "all")
        flnms <- names(fleets)
    Dim <- dim(fleets[[1]]@effort[, years, ])[c(2, 4, 6)]
    Dimnm <- dimnames(fleets[[1]]@effort[, years, ])
    n <- prod(Dim) * length(flnms)
    res <- data.frame(year = rep(years, prod(Dim[2:3]) * length(flnms)),
        quarter = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3] *
            length(flnms)), fleet = rep(flnms, each = prod(Dim)),
        iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(flnms)),
        capacity = numeric(n), costs = numeric(n), effort = numeric(n),
        profits = numeric(n))
    k <- 1
    for (f in flnms) {
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        res[k:(k + prod(Dim) - 1), "capacity"] <- c(fl@capacity[,
            years, ])
        res[k:(k + prod(Dim) - 1), "effort"] <- c(fl@effort[,
            years, ])
        if(!is.null(covars)) res[k:(k + prod(Dim) - 1), "costs"] <- c(costs_flbeia(fl, covars, f)[,years, ])
        res[k:(k + prod(Dim) - 1), "profits"] <- c(revenue_flbeia(fl)[,
            years, ]) - res[k:(k + prod(Dim) - 1), "costs"]
        k <- k + prod(Dim)
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
     return(fleet@fcost*covars[['NumbVessels']][flnm,])            
}


#------------------------------------------------------------------------------#
# catchFlSum data.frame[year, quarter, stock, fleet, iter, ||,|| 
#        landings, discards, price, tacshare] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia
catchFlSum <- function(fleets, advice, flnms = 'all', stknms = 'all', years){
    
    if(flnms[1] == 'all') flnms <- names(fleets)
    if(stknms[1] == 'all') stknms <- catchNames(fleets)
     
    Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years,])
    
    res <- NULL
    
    
                                   
    for(f in flnms){
        
        fl   <- fleets[[f]]

        stfl <- catchNames(fl)        
        sts   <- stknms[stknms %in% stfl]
        
        n <- prod(Dim)*length(sts)
        
        dff <- data.frame(year = rep(years, prod(Dim[2:3])*length(sts)), 
                    quarter = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(sts)), 
                    fleet = rep(f, n), 
                    stock = rep(sts, each = prod(Dim)),
                    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(sts)),  
                    landings = numeric(n), 
                    discards = numeric(n),
                    price    = numeric(n),
                    tacshare = numeric(n))
        
        k <- 1
        
        for(st in sts){
            
            dff[k:(prod(Dim) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,4,6), sum)[,years])    
            dff[k:(prod(Dim) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,4,6), sum)[,years]) 
            dff[k:(prod(Dim) + k-1),'price']    <- c(price_flbeia(fl, st)[,years])
            dff[k:(prod(Dim) + k-1),'tacshare'] <-c((advice$TAC[st,]*advice$quota.share[[st]][f,])[,years])
            
            k <- k + prod(Dim)     
        }
        res <- rbind(res, dff)
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
# catchMtSum data.frame[year, quarter, stock, fleet, metier, iter ||,|| 
#        landings, discards, price] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia 
catchMtSum <- function(fleets, flnms = 'all', stknms = 'all', years){
    
    if(flnms[1] == 'all') flnms <- names(fleets)
    if(stknms[1] == 'all') stknms <- catchNames(fleets)
     
    Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years,])

    res <- NULL
    
    for(f in flnms){
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        for(m in mts){
            mt   <- fl@metiers[[m]]
            stmt <- catchNames(mt)        
            sts  <- stknms[stknms %in% stmt]

            n <- prod(Dim)*length(sts)
        
            dfm <-  data.frame(year = rep(years, prod(Dim[2:3])*length(sts)), 
                        quarter = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(sts)), 
                        fleet = rep(f, n), 
                        metier = rep(m, n),
                        stock = rep(sts, each = prod(Dim)),
                        iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(sts)),  
                        landings = numeric(n), 
                        discards = numeric(n),
                        price = numeric(n))
            k <- 1
            
            for(ss in sts){
                cc <- mt@catches[[ss]]
                dfm[k:(k+prod(Dim)-1),'landings'] <- c(apply(cc@landings[,years,], c(2,4,6), sum, na.rm=T))
                dfm[k:(k+prod(Dim)-1),'discards'] <- c(apply(cc@discards[,years,], c(2,4,6), sum, na.rm=T))
                revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,4,6), sum, na.rm=T)[,years,]
                dfm[k:(k+prod(Dim)-1),'price']  <- c(revst)/dfm[k:(k+prod(Dim)-1),'landings']  
                k <- k + prod(Dim)
            }
            res <- rbind(res, dfm) 
            
        }  
        
                   
    }
    return(res)
}

#------------------------------------------------------------------------------#
# effortMtSum data.frame[year, quarter, fleet, metier, iter ||,|| 
#        effort, effshare] 
#------------------------------------------------------------------------------#
#' @rdname summary_flbeia
effortMtSum <- function(fleets, flnms, years){
    
    if(flnms[1] == 'all') flnms <- names(fleets)
     
    Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years,])

    res <- NULL
    
    for(f in flnms){
        fl <- fleets[[f]]
        mts <- names(fl@metiers)
        n <- prod(Dim)*length(mts)
        
        dff <-  data.frame(year = rep(years, prod(Dim[2:3])*length(mts)), 
                    quarter = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(mts)), 
                    fleet = rep(f, n), 
                    metier = rep(mts, each = prod(Dim)),
                    iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(mts)),  
                    effshare = numeric(n), 
                    effort = numeric(n))
        k <- 1
        for(m in mts){
            mt <- fl@metiers[[m]]
            dff[k:(k+prod(Dim)-1),'effort']   <- c((fl@effort*mt@effshare)[,years,])
            dff[k:(k+prod(Dim)-1),'effshare'] <- c(mt@effshare[,years,])
            k <- k + prod(Dim)
        }
        
        res <- rbind(res, dff)            
    }
    return(res)
}
     
     
           