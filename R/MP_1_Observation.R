#-------------------------------------------------------------------------------
#                           Observation models
#   - observation.mp 
#
# Dorleta GarcYYYa
# Created: 21/12/2010 14:34:47
# Changed: 27/04/2011 09:23:01
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# assessment.mp(stocks, fleets.obs, indices, assess.ctrl)   
#  For the time being only single stock assessments.
#-------------------------------------------------------------------------------
observation.mp <- function(biols, fleets, covars, indices, advice, obs.ctrl, year){

    stnms <- names(biols)
    
    stocks         <- vector('list', length(stnms)) 
    indices.upd    <- indices #vector('list', length(stnms))
    names(stocks) <- stnms
    
# 1 stock per biol
# 1 FLIndices per stock
# 1 FLFleets from fleets the number of fleets/metiers do not need to be equal, 
# for the time being not fleets observation implemented.
     
    for(st in stnms){
        
        # Generete the FLStock
        stocks[[st]] <- eval(call(obs.ctrl[[st]][['stkObs']][['stkObs.model']], biol = biols[[st]], fleets = fleets,  advice = advice,
                        covars = covars, obs.ctrl = obs.ctrl[[st]][['stkObs']], year = year, stknm = st))
        
        # Update the FLIndices [year-1]
        #indices observation (stock by stock)
        nind     <- length(indices[[st]])
        indSt    <- indices[[st]]
        indStnms <- names(indSt)
       # indices.upd[[st]] <- FLIndices()
        
        # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol. 
      #  yrnm.1 <- dimnames(biols[[1]]@n)[[2]][year-1] 
                                                                                               
        for(id in indStnms){
       #     indX <- trim(indices[[st]][[id]], year = dimnames(indices[[st]][[id]]@index)[[2]][1]:yrnm.1)
            obs.model <- obs.ctrl[[st]][['indObs']][[id]][['indObs.model']]
            obs.model <- ifelse(obs.model == 'NoObservation','NoObsIndex', obs.model)
            indices.upd[[st]][[id]] <- eval(call(obs.model, biol = biols[[st]], index = indices[[st]][[id]],
                                            fleets = fleets, covars = covars, obs.ctrl = obs.ctrl[[st]][['indObs']][[id]], year = year))
            
        }
        
        names(indices.upd[[st]]) <- indStnms 
        
    }

    # TO DO: 
    # fleet observaction.
    fleets.obs <- NULL    
    
    return(list(stocks = stocks, fleets.obs = fleets.obs, indices = indices.upd))
}  


NoObsStock <- function(...) return(NULL)

NoObsIndex <- function(...) return(FLIndex())
