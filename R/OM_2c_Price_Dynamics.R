#-------------------------------------------------------------------------------
#                     PRICE DYNAMIC FUNCTIONS
#   Functions to calculate and update the price
#
#    - 'fixedPrice'   - The Price is given as input, just return the object as it is.
#    - 'elasticPrice.total' - elasticity function to model the price (total landings). (Kraak. 2004)
#    - 'elasticPrice.fleet' - elasticity function to model the price (fleets' landings). (Kraak. 2004)
# Dorleta GarcYYYa
# Created: 09/11/2010 09:56:11
# Changed:30/11/2010 12:22:06
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# fixedPrice(fleets, fleets.ctrl, year = 1, season = 1)
#   - Two options: 
#-------------------------------------------------------------------------------
fixedPrice <- function(fleets, covars, fleets.ctrl, year = 1, season = 1,...){
    return(fleets)
}

#-------------------------------------------------------------------------------
# elasticPrice(fleets, covars, fleets.ctrl, stnm, flnm, year = 1, season = 1)
# The price varies according to Kraak elasticity function and are conditioned 
# to the total or fleet's landings. 
# By means of 'covar' we could analyse the effect of imports in a more 
# sophisticated function
#-------------------------------------------------------------------------------

elasticPrice <- function(fleets, covars, fleets.ctrl, stnm, flnm, mtnm, year = 1, season = 1){

    # Parameters
    elas     <- fleets.ctrl[[flnm]][[mtnm]][[stnm]][['pd.els']][,season,] # [na,it] the price and its parameters depend on season.
    La0      <- fleets.ctrl[[flnm]][[mtnm]][[stnm]][['pd.La0']][,season,] # [na,it] the price and its parameters depend on season.
    Pa0      <- fleets.ctrl[[flnm]][[mtnm]][[stnm]][['pd.Pa0']][,season,] # [na,it] the price and its parameters depend on season.
    total    <- fleets.ctrl[[flnm]][[mtnm]][[stnm]][['pd.total']]         # Logic: The function depends on total landings or fleet's landings
    
    fms <- fleets[[flnm]][[mtnm]][[stnm]]
    yr <- year
    ss <- season
    
    # Landings.
    
    if(total == TRUE){
        Lau <- landWStock(fleets, stnm)[,yr,,ss]
#        print('TOTAL')
    }else{
        Lau <- fms@landings.wt 
 #       print('FLEET')
 }
        
    La  <- unitSums(Lau)[drop=T]    # [na,it]
    nu  <- dim(Lau)[3]

    La[La==0] <- 0.001
      
    Pa <- Pa0*(La0/La)^elas    #  [na,it]
    
    # When La = 0 -> Pa = Inf -> set Pa = NA
    Pa <- ifelse( Pa==Inf, NA, Pa)
                                   
    fms@price[,yr,i,ss] <- Pa
    

    fleets[[flnm]]@metiers[[mtnm]]@catches[[stnm]] <- fms
    
    return(fleets)
}
    
