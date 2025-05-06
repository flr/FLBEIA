#-------------------------------------------------------------------------------
#                     PRICE DYNAMIC FUNCTIONS
#            Functions to calculate and update the price
#
#    - 'fixedPrice'      - The Price is given as input, just return the object as it is.
#    - 'elasticPrice'    - Elasticity function to model the price using the same elastic price function for all ages.
#    - 'elasticPriceAge' - Elasticity function to model the price using different elastic price function for each ages.
#
#
### In 'elasticPrice' and 'elasticPriceAge'
#
## Using the 'total' argument, you can choose to condition the price function based on the total landings (TRUE) or on the landings of the fleet in question (FALSE). 
#
## Using the 'type' argument you can choose the elasticity function
#    - 'type1' - The price varies according to Kraak’s elasticity function (Kraak. 2004).
#    - 'type2' - The price varies according to  Isabella’s function 3.
#    - 'type3' - The price varies according to  Isabella’s function 1.
#    - 'type4' - The price varies according to  Isabella’s function 4.
#    - 'type5' - 'fixedPrice'. This option is only available in 'elasticPriceAge'.
#
## Using the 'RP' argument, you can choose to increase or decrease the price by a specific percentage once the new price is calculated. 
#
#
# Dorleta García
# Created: 09/11/2010 09:56:11
# Changed: 30/11/2010 12:22:06
#
# New types of elasticity functions were added by Miren Altuna-Etxabe: 05/07/2024
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# fixedPrice(fleets, covars, fleets.ctrl, year = 1, season = 1)
#-------------------------------------------------------------------------------
fixedPrice <- function(fleets, covars, fleets.ctrl, year = 1, season = 1,...){
    return(fleets)
}

#-------------------------------------------------------------------------------
# elasticPrice(fleets, covars, fleets.ctrl, stnm, flnm, year = 1, season = 1)
#-------------------------------------------------------------------------------

elasticPrice <- function(fleets, covars, fleets.ctrl, stnm, flnm, mtnm, year = 1, season = 1){

    fms      <- fleets[[flnm]][[mtnm]][[stnm]]
    fms.ctrl <- fleets.ctrl[[flnm]][[mtnm]][[stnm]]
    
    yr <- year
    ss <- season
    
    # Parameters 
    elsa     <- fms.ctrl[['pd.elsa']][,ss,,,]                                     # [na,it] Price-landing elasticity.
    La0      <- fms.ctrl[['pd.La0']][,ss,,,]                                      # [na,it] Base landings.
    Pa0      <- fms.ctrl[['pd.Pa0']][,ss,,,]                                      # [na,it] Base price.
    type     <- fms.ctrl[["type"]]                                                # [n,it]  Numeric vector 
    total    <- fms.ctrl[['pd.total']]                                            # Logic:  If TRUE, the function depends on total landings, and if FALSE on the landings of the fleet in question.
    
    if(is.null(fms.ctrl[['pd.RP']])){RP   <- FALSE                                # By default the price function does not have an added increase/decrease.
    } else{  
      RP       <- fms.ctrl[['pd.RP']]                                             # Logic: If TRUE, the function has an added increase/decrease, and if FALSE it does not.
      addRP    <- fms.ctrl[['pd.addRP']]                                          # [n,it] Percentage of the price increased/decreased. Same for all ages.
      yr.addRP <- fms.ctrl[['pd.yr']]                                             # [n,it] Numeric vector. The year from which the price starts increasing or decreasing annually.
    }
    
    # Landings
    if(total == TRUE){
      Lau      <- landWStock(fleets, stnm)[,yr,,ss]
      Lau_init <- landWStock(fleets, stnm)[,yr-1,,ss] 
      Lau_hist <- landWStock(fleets, stnm)[,,,ss] 
    }else{
      Lau      <- fms@landings.wt[,yr,,ss] * fms@landings.n[,yr,,ss]
      Lau_init <- fms@landings.wt[,yr-1,,ss] * fms@landings.n[,yr-1,,ss]
      Lau_hist <- fms@landings.wt[,,,ss] * fms@landings.n[,,,ss]
    }
    
    La_f     <- unitSums(Lau)[drop=T]
    La_init  <- unitSums(Lau_init)[drop=T]
    La_hist  <- unitSums(Lau_hist)[drop=T]
    
    # Type
    if(type == 1) { P <- Pa0*(La0/La_f)^elsa }
    
    if(type == 2 | type == 3 | type == 4) {  
      els    <- quantMeans(elsa)                                                  # [n,it]
      L_f    <- sum(La_f)                                                         # [n,it]
      L_init <- sum(La_init)                                                      # [n,it]
      L_hist <- colSums(La_hist)                                                  # [n,it]
      P_init <- quantMeans(fms@price[, yr-1, , ss])                               # [n,it]
      
      if(L_f == 0){P <- NA                                                        # When L_f = 0 -> set P = NA
      } else{
        if(L_init== 0) {
          L_hist[is.na(L_hist)] <- 0
          filtered_L <- L_hist[L_hist != 0]                                       # Remove years with no landing data.
          Last_Lyr <- max(as.numeric(names(filtered_L)))                          # Identify the last year with landing != 0.
          L_init <- as.numeric(L_hist[names(L_hist)== Last_Lyr])                  # Take the last landing != 0.
          P_init <- quantMeans(fms@price[, grepl(Last_Lyr, names(L_hist)), , ss]) # Take the corresponding price data.
        } 
        if(type == 2){P <- P_init*(L_f/L_init)^els}
        if(type == 3){P <- P_init*(1+els*((L_f-L_init)/L_init))}
        if(type == 4){P <- P_init*exp(els*L_f)}
      }
    }
    
    P <- ifelse(P == Inf, NA, P)                                                  # When L_f = 0 -> P = Inf -> set P = NA
    
    # Added increase/decrease in Price.
    if(RP == TRUE){
      if(type == 1){
        yr.addRP.pos <- which(colnames(La_hist) %in% yr.addRP)                    # Identify the position of the year from which the price starts increasing or decreasing annually.
        nyr <- yr-yr.addRP.pos+1                                                  # Number of years since the first year of the simulation.
      } else{
        P_data <- quantMeans(fms@price[, 1:(yr-1), , ss]) 
        P_data[is.na(P_data)] <- 0
        filtered_P <- P_data[P_data != 0]                                         # Remove years with no price data.
        Last_Pyr <- max(as.numeric(colnames(filtered_P)))                         # Identify the last year with price data != 0.
        
        yr.name <- as.numeric(colnames(La_hist)[yr])                              # Actual year
        nyr <- length(yr.name:Last_Pyr)-1                                         # Number of years since the last price data.
      }
      P <- P*(1+addRP)^nyr
    }
    
    fms@price[, yr, , ss] <- P 
    
    fleets[[flnm]]@metiers[[mtnm]]@catches[[stnm]] <- fms
    
    return(fleets)
}

#-------------------------------------------------------------------------------
# elasticPriceAge(fleets, covars, fleets.ctrl, stnm, flnm, mtnm, year = 1, season = 1)
#-------------------------------------------------------------------------------

elasticPriceAge <- function(fleets, covars, fleets.ctrl, stnm, flnm, mtnm, year = 1, season = 1){
  
    fms      <- fleets[[flnm]][[mtnm]][[stnm]]
    fms.ctrl <- fleets.ctrl[[flnm]][[mtnm]][[stnm]]
    
    yr <- year
    ss <- season
    
    # Parameters
    elsa    <- fms.ctrl[['pd.elsa']][,ss,,,]                                      # [na,it] Price-landing elasticity.
    La0     <- fms.ctrl[['pd.La0']][,ss,,,]                                       # [na,it] Base landings.
    Pa0     <- fms.ctrl[['pd.Pa0']][,ss,,,]                                       # [na,it] Base price.
    type    <- fms.ctrl[["type"]]                                                 # [na,it] Numeric vector: The type depend on age.
    total   <- fms.ctrl[['pd.total']]                                             # Logic: If TRUE, the function depends on total landings, and if FALSE on the landings of the fleet in question.
    
    if(is.null(fms.ctrl[['pd.RP']]))    {RP   <- FALSE                            # By default the price function does not have an added increase/decrease.
    } else{  
      RP       <- fms.ctrl[['pd.RP']]                                             # Logic: If TRUE, the function has an added increase/decrease, and if FALSE it does not.
      addRP    <- fms.ctrl[['pd.addRP']]                                          # [n,it] Percentage of the price increased/decreased. Same for all ages.
      yr.addRP <- fms.ctrl[['pd.yr']]                                             # [n,it] Numeric vector. The year from which the price starts increasing or decreasing annually.
    }
    
    # Landings
    if(total == TRUE){
      Lau      <- landWStock(fleets, stnm)[,yr,,ss]
      Lau_init <- landWStock(fleets, stnm)[,yr-1,,ss] 
      Lau_hist <- landWStock(fleets, stnm)[,,,ss] 
    }else{
      Lau      <- fms@landings.wt[,yr,,ss] * fms@landings.n[,yr,,ss]
      Lau_init <- fms@landings.wt[,yr-1,,ss] * fms@landings.n[,yr-1,,ss]
      Lau_hist <- fms@landings.wt[,,,ss] * fms@landings.n[,,,ss]
    }
    
    La_f    <- unitSums(Lau)[drop=T]
    La_init <- unitSums(Lau_init)[drop=T]
    La_hist <- unitSums(Lau_hist)[drop=T]
    
    Pa_init <- fms@price[, yr-1, , ss]                                            # [na,it]
    
    # Type
    for(a in 1:length(fms@range[["min"]]:fms@range[["max"]])){
      
      tpa <- type[a]                                                              # Price function of age a
      
      if(tpa == 1){Pa <- Pa0[a]*(La0[a]/La_f[a])^elsa[a]
      } else {
        if(La_f[a] == 0){Pa <- NA                                                 # When L_f = 0 -> set P = NA
        } else{
          if(La_init[a]== 0) {
            La_hist[a,][is.na(La_hist[a,])] <- 0
            filtered_La <- La_hist[a,][La_hist[a,] != 0]                          # Remove years with no landing data.
            Last_Layr <- max(as.numeric(names(filtered_La)))                      # Identify the last year with landing != 0.
            La_init <- as.numeric(La_hist[,colnames(La_hist)== Last_Layr])        # Take the last landing =! 0.
            Pa_init <- fms@price[, grepl(Last_Layr, colnames(La_hist)), , ss]     # Take the corresponding price data.
          }
        }
        if(tpa == 2){Pa <- Pa_init[a]*(La_f[a]/La_init[a])^elsa[a]} 
        if(tpa == 3){Pa <- Pa_init[a]*(1+elsa[a]*((La_f[a]-La_init[a])/La_init[a]))} 
        if(tpa == 4){Pa <- Pa_init[a]*exp(elsa[a]*La_f[a])} 
        if(tpa == 5){Pa <- elsa[a]}
      }
      
      Pa <- ifelse( Pa==Inf, NA, Pa)                                              # When La_f = 0 -> Pa = Inf -> set Pa = NA
      
      # Added increase/decrease in Price.
      if(RP == TRUE){
        if(tpa == 1 | tpa == 5){
          yr.addRP.pos <- which(colnames(La_hist) %in% yr.addRP)                  # Identify the position of the year from which the price starts increasing or decreasing annually.
          nyr <- yr-yr.addRP.pos+1                                                # Number of years since the first year of the simulation.
        } else{
          P_data <- quantMeans(fms@price[, 1:(yr-1), , ss]) 
          P_data[is.na(P_data)] <- 0
          filtered_P <- P_data[P_data != 0]                                       # Remove years with no price data.
          Last_Pyr <- max(as.numeric(colnames(filtered_P)))                       # Identify the last year with price data != 0.
          
          yr.name <- as.numeric(colnames(La_hist)[yr])                            # Actual year.
          nyr <- length(yr.name:Last_Pyr)-1                                       # Number of years since the last price data.
        }
        Pa <- Pa*(1+addRP)^nyr
      }
      fms@price[a,yr,,ss] <- Pa                                                   # [na,it]
    }
    
    fleets[[flnm]]@metiers[[mtnm]]@catches[[stnm]] <- fms
    
    return(fleets)
}
