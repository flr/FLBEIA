#------------------------------------------------------------------------------#
#    Auxiliary functions to summarize the results
#
# F_beia(obj)
# SSB_beia(obj)
# B_beia(obj)
# R_beia(obj)
# C_beia(obj);  L_beia(obj); D_beia(obj)
# summary_beia(obj)
#  obj = BEIA output.
#  
#   
# Dorleta García
# Created: 30/01/2011 20:50:27 
# Changed: 30/01/2011 20:50:32
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# F_beia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
F_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){
        # harvest: * if age structured calculate it from 'n'.
        #          * if biomass dyn => assume C = q*E*B => C = F*B and F = C/B.
        na <- dim(obj$biols[[stk]]@n)[1]
        
        if(na == 1){
            # Catch:
            catch <- apply(catchStock(obj$fleets, stk),c(2,6), sum)[drop = TRUE] # [ny,it]
            B     <- (obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,,,1,drop= TRUE] # [ny, it] , 1st season biomass
            res[stk,,] <- catch/B
        }
        else{ 
            fbar_age <- ac(obj$biols[[stk]]@range[c('minfbar')]:obj$biols[[stk]]@range[c('maxfbar')])
            
            Dnms <- list(age = fbar_age, year = 1:ny, iter = 1:it)
            aux  <- array(dim = c(length(fbar_age), ny,it), dimnames = Dnms)           
            
            n.  <- array(unitSums(obj$biols[[stk]]@n)[fbar_age,,,1,drop=T], dim = c(length(fbar_age),ny,it), dimnames = Dnms)                    
            m.  <- array(seasonSums(unitMeans(obj$biols[[stk]]@m))[fbar_age,drop=T], dim = c(length(fbar_age),ny,it), dimnames = Dnms)           
            c.  <- array(apply(catchStock(obj$fleets, stk),c(1:2,6), sum)[fbar_age,drop = TRUE], dim = c(length(fbar_age),ny,it), dimnames = Dnms) 
        
            fobj <- function(f,n,m,c){ return( f/(f+m)* (1-exp(-(f+m)))*n -c)}
        
            for(y in 1:ny){
                for(a in fbar_age){
                    for(i in 1:it){
                        aux[a,y,i] <- ifelse(n.[a,y,i] == 0, 0,
                                                uniroot(fobj, lower = 0, upper = 10, n = n.[a,y,i], m=m.[a,y,i], c = c.[a,y,i])$root)  
            }}}
           res[stk,,] <- apply(aux,2:3,mean) 
        }
    }
    return(res)
}


#------------------------------------------------------------------------------#
# SSB_beia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
SSB_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # SSB 1st season
        res[stk,,] <- apply(unitSums(obj$biols[[stk]]@n*obj$biols[[stk]]@wt*obj$biols[[stk]]@fec)[,,,1], c(2,6), sum, na.rm = TRUE)[drop=T]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# B_beia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
B_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(unitSums(obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,,,1], c(2,6), sum, na.rm = TRUE)[drop=T]
    }
    return(res)
}

#------------------------------------------------------------------------------#
# R_beia(obj) :: res[stocks, years, it] 
# If age struc => recruitment.
# If biodyn    => growth.
#------------------------------------------------------------------------------#
R_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # 
        na <- dim(obj$biols[[stk]]@n)[1]
        if(na > 1){
            res[stk,,] <- obj$biols[[stk]]@n[1,,1,1,drop=T]
            if(dim(obj$biols[[stk]]@n)[3]>1){
                for(ss in 2:dim(obj$biols[[stk]]@n)[3]) res[stk,,] <- res[stk,,] + obj$biols[[stk]]@n[1,,ss,ss,drop=T]
            }
        }else{
            catch <- matrix(apply(catchStock(obj$fleets, stk),c(2,6), sum)[drop = TRUE],ny,it) # [ny,it]
            B     <- matrix((obj$biols[[stk]]@n*obj$biols[[stk]]@wt)[,,,1,drop= TRUE],ny,it) # [ny, it] , 1st season biomass
            res[stk,-ny,] <- B[-1,] - B[-ny,] + catch[-ny,]
            
        }
    }
    return(res)
}

#------------------------------------------------------------------------------#
# C_beia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
C_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(catchWStock(obj$fleets, stk),c(2,6), sum)[drop = TRUE] # [ny,it]
    }
    return(res)
}

#------------------------------------------------------------------------------#
# L_beia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
L_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(landWStock(obj$fleets, stk),c(2,6), sum)[drop = TRUE] # [ny,it]
    }
    return(res)
}

#------------------------------------------------------------------------------#
# D_beia(obj) :: res[stocks, years, it] 
#------------------------------------------------------------------------------#
D_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it), dimnames = list(stock = stknms, year = yrnms))
    
    for(stk in stknms){ # B 1st season
        res[stk,,] <- apply(discWStock(obj$fleets, stk),c(2,6), sum)[drop = TRUE] # [ny,it]
    }
    return(res)
}


#------------------------------------------------------------------------------#
# summary_beia(obj) :: res[stocks, years, it, indicators] 
#------------------------------------------------------------------------------#
summary_beia <- function(obj){
    stknms <- names(obj$biols)
    
    it     <- dim(obj$biols[[1]]@n)[6]
    ny     <- dim(obj$biols[[1]]@n)[2]
    yrnms  <- dimnames(obj$biols[[1]]@n)[[2]]
    
    res <- array(dim = c(length(stknms), ny,it, 7), dimnames = list(stock = stknms, year = yrnms, iter = 1:it, 
                                                      indicators = c('rec', 'ssb', 'f', 'biomass', 'catch', 'landings', 'discards')))
    
    res[,,,1] <- R_beia(obj)
    res[,,,2] <- SSB_beia(obj)
    res[,,,3] <- F_beia(obj)
    res[,,,4] <- B_beia(obj)
    res[,,,5] <- C_beia(obj)
    res[,,,6] <- L_beia(obj)
    res[,,,7] <- D_beia(obj)

    return(res)
    
}


#------------------------------------------------------------------------------#
# BIOsummary(obj) :: DATA.FRAME[stocks, years, it, indicators, value] 
#------------------------------------------------------------------------------#
bioSum <- function(object){
    xx <- summary_beia(object)
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
ecoSum <- function(fleets, flnms = 'all', years){
    
    if(flnms[1] == 'all') flnms <- names(fleets)
    
    Dim   <- dim(fleets[[1]]@effort[,years,])[c(2,4,6)]
    Dimnm <- dimnames(fleets[[1]]@effort[,years,])
    
    n    <- prod(Dim)*length(flnms)
    
    res <- data.frame(year = rep(years, prod(Dim[2:3])*length(flnms)), 
                      quarter = rep(rep(Dimnm[[4]], each = Dim[1]), Dim[3]*length(flnms)), 
                      fleet = rep(flnms, each = prod(Dim)), 
                      iter = rep(rep(1:Dim[3], each = prod(Dim[1:2])), length(flnms)),  
                      capacity = numeric(n), 
                      costs = numeric(n), 
                      effort = numeric(n),
                      profits = numeric(n))
    k <- 1  
                           
    for(f in flnms){
        
        fl   <- fleets[[f]]
        mts  <- names(fl@metiers)
        
        res[k:(k+prod(Dim)-1),'capacity'] <- c(fl@capacity[,years,])
        res[k:(k+prod(Dim)-1),'effort']   <- c(fl@effort[,years,])
        res[k:(k+prod(Dim)-1),'costs']    <- c(costs_beia(fl)[,years,])
        res[k:(k+prod(Dim)-1),'profits']  <- c(revenue_beia(fl)[,years,]) -  res[k:(k+prod(Dim)-1),'costs']
        
        k <- k + prod(Dim)
    }
    return(res)
}                               

    


#-------------------------------------------------------------------------------
# revenue_beia(fleet, years)
#-------------------------------------------------------------------------------
revenue_beia <- function(fleet){
    
    sts <- catchNames(fleet)
    mts <- names(fleet@metiers)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        m <- fleet@metiers[[mt]]
        for(st in sts){
            if(!(st %in% catchNames(m))) next
            dat <- m@catches[[st]]
            res <- res + apply(dat@landings.n*dat@landings.wt*dat@price, c(2,4,6),sum)
        }
    }
    return(res)               
}

#-------------------------------------------------------------------------------
# costs_beia(fleet, years)
#-------------------------------------------------------------------------------
costs_beia <- function(fleet){
    
    res <- totvcost_beia(fleet) + totfcost_beia(fleet)
    
    return(res)               
}

#-------------------------------------------------------------------------------
# totvcost_beia(fleet, years)
#-------------------------------------------------------------------------------
totvcost_beia <- function(fleet){
    
    mts <- names(fleet@metiers)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        res <- res + fleet@metiers[[mt]]@vcost*fleet@effort*fleet@metiers[[mt]]@effshare
    }
    return(res)               
}

#-------------------------------------------------------------------------------
# totvcost_beia(fleet, years)
#-------------------------------------------------------------------------------
totfcost_beia <- function(fleet){
     return(fleet@fcost*fleet@capacity)            
}


#------------------------------------------------------------------------------#
# catchFlSum data.frame[year, quarter, stock, fleet, iter, ||,|| 
#        landings, discards,price] 
#------------------------------------------------------------------------------#
catchFlSum <- function(fleets, flnms = 'all', stknms, years){
    
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
                    price    = numeric(n))
        
        k <- 1
        
        for(st in sts){
            
            dff[k:(prod(Dim) + k-1),'landings'] <- c(apply(landWStock.f(fl, st),c(2,4,6), sum)[,years])    
            dff[k:(prod(Dim) + k-1),'discards'] <- c(apply(discWStock.f(fl, st),c(2,4,6), sum)[,years]) 
            dff[k:(prod(Dim) + k-1),'price']    <- c(price_beia(fl, st)[,years])
            
            k <- k + prod(Dim)     
        }
        res <- rbind(res, dff)
    }
    return(res)
}                               


#-------------------------------------------------------------------------------
# price_beia(fleet, years)(mean price in a fleet)
#-------------------------------------------------------------------------------
price_beia <- function(fleet, stock){

    mts <- names(fleet@metiers)
    
    totL <- apply(landWStock.f(fleet, stock), c(2,4,6), sum)
    
    res <- FLQuant(0, dimnames = dimnames(fleet@effort))
    
    for(mt in mts){
        m <- fleet@metiers[[mt]]
        if(!(stock %in% catchNames(m))) next
        dat <- m@catches[[stock]]
        res <- res + apply(dat@landings.n*dat@landings.wt*dat@price, c(2,4,6),sum)
    }
    
    res <- res/totL
    
    return(res)                
}

#------------------------------------------------------------------------------#
# effortMtSum data.frame[year, quarter, fleet, metier, iter ||,|| 
#        effort, effshare] 
#------------------------------------------------------------------------------#
 
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
                dfm[k:(k+prod(Dim)-1),'landings'] <- c(apply(cc@landings[,years,], c(2,4,6), sum))
                dfm[k:(k+prod(Dim)-1),'discards'] <- c(apply(cc@discards[,years,], c(2,4,6), sum))
                revst <- apply(cc@landings.n*cc@landings.wt*cc@price, c(2,4,6), sum)[,years,]
                dfm[k:(k+prod(Dim)-1),'price']  <- c(revst)/dfm[k:(k+prod(Dim)-1),'landings']  
                k <- k + prod(Dim)
            }
            res <- rbind(res, dfm) 
            
        }  
        
                   
    }
    return(res)
}



#------------------------------------------------------------------------------#
# catchMtSum data.frame[year, quarter, stock, fleet, metier, iter ||,|| 
#        landings, discards, price] 
#------------------------------------------------------------------------------#

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
     
     
           