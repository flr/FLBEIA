#-------------------------------------------------------------------------------
#                   SCD :: Simple Capital Dynamics 
#
# Dorleta Garcia
# created: 04/06/2012 15:52:25
# changed: 18/06/2012 10:20:52
#-------------------------------------------------------------------------------

# Crewcost: FixedCost(Salaries) + VariableCost(CrewShare).
# Crewshare = % of the total landing value that belongs to the crew.

SCD <- function(fleets, covars, fleets.ctrl, flnm, year = 1, season = 1,...){
    
   
    
    fleet <- fleets[[flnm]]
    
    ny <- dim(fleet@effort)[2]
    ns <- dim(fleet@effort)[4]
    it <- dim(fleet@effort)[6]
    
    # VaC
    VaC <- seasonSums(totvcost_flbeia(fleet)[,year]) # total anual variable costs
    # FxC
    FxC <- seasonSums(covars[["NumbVessels"]][flnm, ] * fleet@fcost)[, year]
    # FuC  # per unit of effort, we asume common cost for all the metiers.
    FuC <- seasonSums(covars[['FuelCost']][flnm,]*fleet@effort)[,year]
    # CaC # per unit of capacity
    CaC <- seasonMeans((covars[['CapitalCost']][flnm,]*covars[["NumbVessels"]][flnm, ]))[,year]
    # Revenue
    Rev <- revenue_flbeia(fleet)[,year]
    Rev <- ifelse(Rev == 0, 1e-16, Rev)
    # CrC
    CrC <- seasonSums((Rev*fleet@crewshare[,year]  +  covars[['Salaries']][flnm,year]))
    
    Rev <- seasonSums(Rev)
    
    x1 <- FuC/Rev
    x2 <- VaC/Rev
    
    a <- CrC + FxC + CaC
    b <- 1 - x1 - x2
    
    BER <- a/b
    
    
    Inv <- c((Rev - BER)/Rev)*c(covars[['InvestShare']][flnm,year,,ns]) # The share in last season
    
    Inv <- ifelse((Rev - BER) < 0 & Rev < 0, -Inv, Inv) # If both are negative, the ratio is positive!! change it!!
    
    Ks <- fleet@capacity[,year][drop=T]    # seasonal capacity [ns,ni]
    K  <- c(seasonSums(fleet@capacity[,year])) # annual capacity. [ni]

    # pKs How annual capacity is distributed along seasons.
    if(ns == 1)      pKs <- rep(1,it) #[ni]
    else  if(it > 1) pKs <- sweep(Ks,2,K,"/")    # ns > 1 [ns,ni]
          else       pKs <- Ks/K    # [ns]
    
    w1 <- c(covars[['w1']][flnm,year,,ns]) # last season value
    w2 <- c(covars[['w2']][flnm,year,,ns]) # last season value 
    
    
#    # Translate Inv in number of vessels.
#    Inv_ves <- ifelse(Inv>0, Inv/c(covars[['NewVessPrice']][flnm, year,]), Inv/c(covars[['OldVessPrice']][flnm, year,]))
  
    omega <- ifelse(Inv < 0, 
                        ifelse(-Inv < w1, Inv*K, -w1*K),  # Inv < 0
                        ifelse(Inv < w2, Inv*K, w2*K))    # Inv >= 0  
                        
  #  print(omega)
                
    # Investment in new vessels only occur if the operational days of existing vessesl is equal to capacity and investment saving is >0.
    # In iters where effort == capacity?    
    # In iterSel although the money for investment is >0 there is no investment.
    Ef <- c(seasonSums(fleet@effort[,year]))
    iterSel <- which(omega > 0 & Ef < 0.99*K)              
    
    omega[iterSel] <- 0
    
    # If year is not last year Update capacity  in year [year+1].
    if (year < ny){
        fleets[[flnm]]@capacity[, year + 1] <- Ks + omega*pKs
        covars[['NumbVessels']][flnm,year+1,] <-   fleets[[flnm]]@capacity[, year + 1]/(covars[['MaxDays']][flnm,year+1,])
    }

    
    return(list(fleets = fleets, covars = covars))
}

