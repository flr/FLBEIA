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
    it <- dim(fleet@effort)[6]
    
    # VaC
    VaC <- totvcost_beia(fleet)[,year]
    # FxC
    FxC <- totfcost_beia(fleet)[,year]
    # FuC  # per unit of effort, we asume common cost for all the metiers.
    FuC <- (covars[['FuelCost']][flnm,]*fleet@effort)[,year]
    # CaC # per unit of capacity
    CaC <- (covars[['CapitalCost']][flnm,]*fleet@capacity)[,year]
    # Revenue
    Rev <- revenue_beia(fleet)[,year]
    # CrC
    CrC <- (Rev*fleet@crewshare[,year])  +  covars[['Salaries']][flnm,year]
    # w1, w2
    w1 <- covars[['w1']][flnm, year,]
    w2 <- covars[['w2']][flnm, year,]
    
    x1 <- FuC/Rev
    x2 <- VaC/Rev
    
    a <- CrC + FxC + CaC
    b <- 1 - x1 - x2
    
    BER <- a/b
    
    
    Inv <- c((Rev - BER)/Rev)*c(covars[['InvestShare']][flnm,year])
    
    K <- c(fleet@capacity[,year]) #capacity.
    
    w1 <- c(covars[['w1']][flnm,year]) 
    w2 <- c(covars[['w2']][flnm,year]) 
    
    
#    # Translate Inv in number of vessels.
#    Inv_ves <- ifelse(Inv>0, Inv/c(covars[['NewVessPrice']][flnm, year,]), Inv/c(covars[['OldVessPrice']][flnm, year,]))
  
    omega <- ifelse(Inv < 0, 
                        ifelse(-Inv < w1, Inv*K, -w1*K),  # Inv < 0
                        ifelse(Inv < w2, Inv*K, w2*K))    # Inv >= 0  
                        
    print(omega)
                
    # Investment in new vessels only occur if the operational days of existing vessesl is equal to capacity and investment saving is >0.
    # In iters where effort == capacity?    
    # In iterSel although the money for investment is >0 there is no investment.
    Ef <- c(fleet@effort[,year])
    iterSel <- which(omega > 0 & Ef < 0.99*K)              
    
    omega[iterSel] <- 0
    
    # If year is not last year Update capacity  in year [year+1].
    if(year < ny) fleets[[flnm]]@capacity[,year +  1] <-  K + omega

    
    return(list(fleets = fleets, covars = covars))
}




