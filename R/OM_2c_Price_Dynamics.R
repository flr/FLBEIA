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

elasticPrice <- function(fleets, covars, fleets.ctrl, stnm, flnm, year = 1, season = 1){

    # Parameters
    elas     <- fleets.ctrl[[flnm]][[stnm]][['pd.els']][,season,] # [na,it] the price and its parameters depend on season.
    La0      <- fleets.ctrl[[flnm]][[stnm]][['pd.La0']][,season,] # [na,it] the price and its parameters depend on season.
    Pa0      <- fleets.ctrl[[flnm]][[stnm]][['pd.Pa0']][,season,] # [na,it] the price and its parameters depend on season.
    total    <- fleets.ctrl[[flnm]][[stnm]][['pd.total']]         # Logic: The function depends on total landings or fleet's landings
    
    f <- fleets[[flnm]]
    yr <- year
    ss <- season
    
    # Landings.
    
    if(total == TRUE){
        Lau <- landWStock(fleets, stnm)[,yr,,ss]
#        print('TOTAL')
    }else{
        Lau <- landWStock.f(f, stnm)[,yr,,ss]   
 #       print('FLEET')
 }
        
    La  <- unitSums(Lau)[drop=T]    # [na,it]
    nu  <- dim(Lau)[3]

    Pa <- Pa0*(La0/La)^elas    #  [na,it]
    
    # When La = 0 -> Pa = Inf -> set Pa = NA
    Pa <- ifelse( Pa==Inf, NA, Pa)
                                   
    for(mt in 1:length(f@metiers)){
        
        if(!(stnm %in% names(f@metiers[[mt]]@catches))) next
        
         for(i in 1:nu)
            f@metiers[[mt]]@catches[[stnm]]@price[,yr,i,ss] <- Pa
    }

    fleets[[flnm]] <- f
    return(fleets)
}
    
