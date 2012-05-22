#-------------------------------------------------------------------------------
#           OBSERVATION MODEL FUNCTIONS
#   - gmeanHistGrowth(stock, y0,y1)
#   - fwdBD(stock,fwdControl)       
#   - biomass.tg, catch.tg, discards.tg, f.tg, fdisc.tg, fland.tg, landings.tg
#
# Dorleta GarcYYYa
# Created: 15/12/2010 19:58:56
# Changed: 16/12/2010 12:35:50
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# fwdBD(stock, ctrl, growth.years)
#    - growth.years: Year range used to calculate population growth to be used in the forecast. 
#-------------------------------------------------------------------------------

fwdBD <- function(stock, ctrl, growth.years)
    {
    
    options(warn = -1)
    if (!validObject(stock))
       stop("stock object not valid")
       
    if (!validObject(ctrl))
       stop("ctrl not valid")

    yrs <- sort(as.numeric(as.character(unique(ctrl@target[, "year"])))) 

    ## check years
    ## years in ctrl have to be in order
    if (!all(yrs==sort(yrs)))
       stop("years in ctrl not in order")
    ## no gaps in years
    if (length(min(yrs):max(yrs))!=length(unique(yrs)))
       stop("years in ctrl not contiguous")
    ## years in ctrl have to be in FLStock stock
    if (!all(ac(yrs) %in% ac(dims(stock)$minyear:(dims(stock)$maxyear))))
       stop("years in ctrl outside of those in stock object")
    ## Need year+1 in FLStock object
    if (max(yrs) == dims(stock)$maxyear){
       endYr<-dims(stock)$maxyear+1
       stock <- window(stock, end=dims(stock)$maxyear+1)}
    else
       endYr<-NULL
    
    ## check iters in ctrl are '1 or n' and correct if necessary
#    ctrl@trgtArray <- chkBDTrgtArrayIters(stock,ctrl@trgtArray,sr)

    ctrl@target    <- chkBDTargetQuantity(ctrl@target,object)

    stock.n(stock)[1,ac(min(ctrl@target[,"year"]))]<-NA

    for(yr in yrs){
        y <- yr - dims(stock)$minyear + 1  # yr position.
         
        control <-   ctrl@target[ctrl@target$year == yr, ]
         
        # There must ONE and ONLY ONE target. 
        if(all(is.na(control$val)))  stop('Target is missing in year: ', yr)
        if(sum(!is.na(control$val)) > 1) stop('Only one target per year is allowed. Several targets in year: ', yr) 
        # In each row a CONSTRAINT OR A TARGET but NOT BOTH.
        mi <- !is.na(control$min)
        ma <- !is.na(control$max)     
        va <- which(!is.na(control$val)) # WE KNOW there is only 1 with TRUE.
        if(any(c(mi[va], ma[va]))) stop('Constraints and targets can not be mixed. Mixed in year : ', yr)
     
        gr <- gmeanHistGrowth(stock, growth.years)
         
        # Fist apply the TARGET, and then the CONSTRAINTS in order of appearance.
        # TARGET's ROW.
        tg     <- control[va,]
        tg.fun <- paste(as.character(tg[1,"quantity"]), "tg", sep = ".")
        tg.val <- as.numeric(tg[1,'val'])
        # If the target is relative calculate the absolute value.   
        if(!is.na(tg['rel.year'])){   # => IT IS RELATIVE. EXTRACT the QUANTITY in THAT YEAR and multitply by the RELATIVE target.
            yrel   <- as.numeric((tg['rel.year'] - dims(stock)$minyear + 1))
            quan   <- as.character(tg[1,'quantity'])
            tg.val <- getQuantity(stock, quan, yrel)*tg.val
        }
        
        # Apply target function
        res.tg <- eval(call(tg.fun, stock = stock , gr = gr, target = tg.val, year = y))
        
        # Update Stock.
        stock@stock.n[,y,]      <- res.tg$biomass/stock@stock.wt[,y,]
        stock@stock[,y,]        <- res.tg$biomass
        stock@stock.n[,y+1,]    <- res.tg$biomass_next/stock@stock.wt[,y,]
        stock@stock[,y+1,]      <- res.tg$biomass_next
        stock@catch.n[,y,]      <- res.tg$catch/stock@catch.wt[,y,]
        stock@catch[,y,]        <- res.tg$catch
        stock@landings.n[,y,]   <- res.tg$landings/stock@landings.wt[,y,]
        stock@landings[,y,]     <- res.tg$landings
        stock@discards.n[,y,]   <- res.tg$discards/stock@discards.wt[,y,]
        stock@discards[,y,]     <- res.tg$discards
        stock@harvest[,y,]      <- res.tg$f
                
        # CONSTRAINTS.
        cnst <- control[-va,]  # data frame with the constraints.
        if(dim(cnst)[1] == 0) next
        
        for(cn in 1:dim(cnst)[1]){
            cnst.qnt <- as.character(cnst[cn,'quantity'])   # c("biomass","catch","landings","discards","f", "f.landings","f.discards") 
            cnst.min <- ifelse(is.na(cnst[cn,'min']), 0, as.numeric(cnst[cn,'min']))    
            cnst.max <- ifelse(is.na(cnst[cn,'max']), Inf, as.numeric(cnst[cn,'max'])) 
            # current value  of 'quantity' (absolute value)
            quan <- getQuantity(stock, cnst.qnt, y)
            # If the constraints are relative we have to extract the 'absolute values'.
            if(!is.na(cnst[cn,'rel.year'])){   # => IT IS RELATIVE. EXTRACT the QUANTITY in THAT YEAR and multitply by the RELATIVE target.
                yrel <- cnst[cn,'rel.year'] - dims(stock)$minyear + 1  
                yy <- getQuantity(stock, cnst.qnt, yrel)
                cnst.min <- yy*cnst.min
                cnst.max <- yy*cnst.max
            }
            if(quan > cnst.min & quan < cnst.max) next
            
            new.tg     <- ifelse(quan < cnst.min, cnst.min, cnst.max)
            new.tg.fun <- paste(cnst.qnt, "tg", sep = ".")
            
            # Apply target function
            res.new.tg <- eval(call(new.tg.fun, stock = stock , gr = gr, target = new.tg, year = y))
        
            # Update Stock.
            stock@stock.n[,y,]      <- res.new.tg$biomass/stock@stock.wt[,y,]
            stock@stock[,y,]        <- res.new.tg$biomass
            stock@stock.n[,y+1,]    <- res.new.tg$biomass_next/stock@stock.wt[,y,]
            stock@stock[,y+1,]      <- res.new.tg$biomass_next
            stock@catch.n[,y,]      <- res.new.tg$catch/stock@catch.wt[,y,]
            stock@catch[,y,]        <- res.new.tg$catch
            stock@landings.n[,y,]   <- res.new.tg$landings/stock@landings.wt[,y,]
            stock@landings[,y,]     <- res.new.tg$landings
            stock@discards.n[,y,]   <- res.new.tg$discards/stock@discards.wt[,y,]
            stock@discards[,y,]     <- res.new.tg$discards
            stock@harvest[,y,]      <- res.new.tg$f
                  
        }                                 
     } 

     options(warn = 0)
    return(stock)
}
    
    
#-------------------------------------------------------------------------------
#  gmeanHistGrowth: Calculate the geometric mean growth to be used for the 
#       projection in fwd.
#-------------------------------------------------------------------------------

gmeanHistGrowth <- function(stock, years){   # y0:y1, pueden ser numericos o caracteres.
    it <- dim(stock@m)[6]
    ny <- length(years)
    
    bio <- matrix(quantSums(stock@stock.n*stock@stock.wt)[,years,drop=T],ny,it)  # matrix[ny,it]
    cw  <- matrix(computeCatch(stock)[,years,drop=T],ny,it)                      # matrix[ny,it]
    
    # The growth is calculated in percemtage to avoid the problem with negative values 
    # in the calculation of geometric  mean.
    growth <- matrix((bio[-1,]+ cw[-ny,])/ bio[-ny,] , ny-1,it)
    
    res <- apply(growth,2, function(x) prod(x)^(1/(ny-1)))
    
    # The percentage is applied to last year biomass.
    res <- bio[ny]*(res - 1)  # res > 1 => positive growth ** res < 1 => negative growth. 
    
    return(res)
    
}


#-------------------------------------------------------------------------------
#  getQuantity: Extract 'quantity' from a FLStockobject
#   quantity = "biomass","catch","landings","discards","f", "f.landings","f.discards" 
#   year = single year ', 'POSITION!!)
#-------------------------------------------------------------------------------

getQuantity <- function(stock, quantity, year){

    if(quantity %in% c("catch","landings","discards")) return(c(slot(stock, quantity)[,year,])) 
    if(quantity == 'biomass')  return(c((stock@stock.wt*stock@stock.n)[,year,]))
    if(quantity == 'f') return(stock@harvest[,year,])
    dr <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T]
    if(quantity == 'f.landings') return(c(stock@harvest*(1-dr)))
    
    return(c(stock@harvest*dr))  # tg['quantity'] == 'f.discards'
  
}  
  
#-------------------------------------------------------------------------------
#  KKK.tg(stock, gr, KKK, year): Project the population based on targets.
#      - stock: FLStock object with na = 1.
#      -    gr: Growth of the population in 'year'
#      -   KKK: The constraint 'numeric(it)'.
#      -  year: Projected year - NUMERIC.
#-------------------------------------------------------------------------------

# When the target is biomass, it is the desired biomass in [year+1] because
# the biomass in [year] is the product of what happens in [year-1].
 
biomass.tg <- function(stock, gr, target, year){ # bio[year+1] -> catch (land, disc)?, f?
    
    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
    Cy   <- By + gr - target                              # [it]
    Cy   <- ifelse(Cy < 0, 0, Cy)
    Fy   <- Cy/target
    
    By1   <- By + gr - Cy
    
    # Discards and landings (3 years average for the ratio)
    dr <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]
    
    return(list(catch = Cy, discards = Cy*dr, landings = Cy*(1-dr), f = Fy, biomass = By, biomass_next = By1))
}

catch.tg <- function(stock, gr, target, year){ #
    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
    By1   <- By + gr - target
    
    catch <- ifelse(target > By, 0.9*By, target)
    Fy <- catch/By                            # [it]
    
    # Discards and landings (3 years average for the ratio)
    dr <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]
    
    return(list(catch = catch, discards = catch*dr, landings = catch*(1-dr), f = Fy, biomass = By, biomass_next = By1))                                     
}



landings.tg <- function(stock, gr, target, year){

    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
    
    # Discards and landings (3 years average for the ratio)
    dr    <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]
    catch <- target/(1-dr)
    
    By1   <- By + gr - catch   
    
    catch <- ifelse(catch > By, 0.9*By, catch)
    Fy <- catch/By                            # [it]

    return(list(catch = catch, discards = catch*dr, landings = catch*(1-dr), f = Fy, biomass = By, biomass_next = By1))      
}

discards.tg <- function(stock, gr, target, year){
    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
    
    # Discards and landings (3 years average for the ratio)
    dr    <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]
    catch <- target/dr
    
    catch <- ifelse(catch > By, 0.9*By, catch)
    Fy <- catch/By                            # [it]
    
    By1   <- By + gr - catch

    return(list(catch = catch, discards = catch*dr, landings = catch*(1-dr), f = Fy, biomass = By, biomass_next = By1))   
}

f.tg <- function(stock, gr, target, year){
    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
    
    Fy <- target
    
    if(Fy >= 1) stop('Harvest rate must be < 1')
    
    catch <- Fy*By
    
    By1   <- By + gr - catch
    
    # Discards and landings (3 years average for the ratio)
    dr    <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]

    return(list(catch = catch, discards = catch*dr, landings = catch*(1-dr), f = Fy, biomass = By, biomass_next = By1))   
}

f.landings.tg <- function(stock, gr, target, year){
    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
   
    # Discards and landings (3 years average for the ratio)
    dr    <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]

    Fy <- target/(1-dr) 
    
    if(Fy >= 1) stop('Harvest rate must be < 1')
    
    catch <- Fy*By
    
    By1   <- By + gr - catch
    
    return(list(catch = catch, discards = catch*dr, landings = catch*(1-dr), f = Fy, biomass = By, biomass_next = By1))   
}   

f.discards.tg <- function(stock, gr, target, year){
    By_1 <- (stock@stock.n*stock@stock.wt)[,year-1,drop=T] # [it]
    Cy_1 <- (stock@catch)[,year-1,drop = T]                # [it]
    By   <- By_1 + gr - Cy_1                               # [it]
   
    # Discards and landings (3 years average for the ratio)
    dr    <- yearMeans((stock@discards/stock@catch)[,(year-3):(year-1),])[drop=T] # [it]

    Fy <- target/dr
    
    if(Fy >= 1) stop('Harvest rate must be < 1')
    
    catch <- Fy*By
    
    By1   <- By + gr - catch
    
    return(list(catch = catch, discards = catch*dr, landings = catch*(1-dr), f = Fy, biomass = By, biomass_next = By1))   
}