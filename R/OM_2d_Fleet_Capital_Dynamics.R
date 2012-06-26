#-------------------------------------------------------------------------------
#                     INVESTMENT DYNAMIC FUNCTIONS
#  * This functions could change capacity and catchability slots.
#  * They work/update only one fleet in each call to the function. 
#
#    - 'fixedCapital' - The capacity and catchability are given as input, 
#                   just return the object as it is.
#
# Dorleta GarcYYYa
# Created: 01/12/2010 11:43:09
# Changed:01/12/2010 11:43:04
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# fixedCapital(fleets, fleets.ctrl, flnm, year = 1, season = 1) 
#-------------------------------------------------------------------------------
fixedCapital <- function(fleets, covars, fleets.ctrl, flnm, year = 1, season = 1,...){
    return(list(fleets = fleets, covars = covars))
}
