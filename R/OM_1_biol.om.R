#-------------------------------------------------------------------------------
#       biols.om(biols, fleets, srs, biols.ctrl, year, season,... )
# 
#   Output: Updated FLBiols 
# 
#   Projects the biological populations from: - [year-1, ns] to [year,1] 
#                                          or, - [year, season-1] to [year,season]
# 
#   Age structured or biomass dynamic populations.
#       * Age structured  => ASP
#       * Biomass dynamic => FLBDsim
#
# Dorleta Garcia
# Created: 25/10/2010 15:43:38
# Changed: 27/10/2010 12:34:33
#-------------------------------------------------------------------------------

biols.om <- function(biols, fleets, SRs, BDs, covars, biols.ctrl, year, season){

    stnms <- names(biols)
    
    for(st in stnms){
        # population dynamic model
        dyn.model <- biols.ctrl[[st]]$dyn.model   
        
        res <- eval(call(dyn.model, biols = biols, fleets = fleets, SRs = SRs, BDs = BDs, stknm = st, year = year, season = season, ctrl = biols.ctrl, covars = covars)) 
         
        biols[[st]] <- res$biol
        
        if(!is.null(SRs[[st]]))  SRs[[st]] <- res$SR
        if(!is.null(BDs[[st]]))  BDs[[st]] <- res$BD
    }
    
    return(list(biols = biols, SRs = SRs, BDs = BDs))

}










