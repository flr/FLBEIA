#-------------------------------------------------------------------------------
#                           Assessments
#   - assessment.mp:
#         update, stock.n, stock and harvest for each stock, 
#   - NoAssessment.
#
# Dorleta GarcYYYa
# Created: 21/12/2010 07:55:33
# Changed: 21/12/2010 07:55:38
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# assessment.mp(stocks, fleets.obs, indices, assess.ctrl)   
#  For the time being only single stock assessments.
#-------------------------------------------------------------------------------
assessment.mp <- function(stocks, fleets.obs, indices, assess.ctrl, datayr){

    stnms <- names(stocks)
     
    for(st in stnms){
        
        if(assess.ctrl[[st]] == "NoAssessment")  next
        
        
        # trim the indices, fron index specific initial yeaer to the assessment year.
        indST <- indices[[st]]
        if(!is.null(indST))
            indST <- FLIndices(lapply(indices[[st]],function(x) trim(x, year = dimnames(x@index)[[2]][1]:datayr)))
        
        if(assess.ctrl[[st]]$work_w_Iter == TRUE){  # The assessment model works with iters => all the iterations are 'adjusted' in one run.
            res <- eval(call(assess.ctrl[[st]][['assess.model']], stock = stocks[[st]], indices = indST, control = assess.ctrl[[st]]$control))
            stock.n(stocks[[st]]) <- res@stock.n
            harvest(stocks[[st]]) <- res@harvest
           # stock(stocks[[st]])   <- res@stock
        }
        else{  # The assessment model _DOES NOT_ work with iters =>  the iterations are 'adjusted' one by one.
            it <- dim(stocks[[st]]@m)[6]
            
            for(i in 1:it){ 
                res <- eval(call(assess.ctrl[[st]][['assess.model']], stock = iter(stocks[[st]],i), 
                            indices = iter(indST,i), control = assess.ctrl[[st]]$control))
                iter(stock.n(stocks[[st]]),i) <- res@stock.n
                iter(harvest(stocks[[st]]),i) <- res@harvest
            }
        }
        stock(stocks[[st]])   <- computeStock(stocks[[st]])
        units(harvest(stocks[[st]])) <- assess.ctrl[[st]]$harvest.units
    }
    
    return(stocks)

}  

#-------------------------------------------------------------------------------
# NoAssessment(stock)   
#-------------------------------------------------------------------------------
NoAssessment <- function(stock,...){
   
    return(stock) 

}  


#-------------------------------------------------------------------------------
#  XSA FOR SIMULATION
#-------------------------------------------------------------------------------
FLXSAnew <- function(stock, indices, control) return(FLXSA(stock, indices, control, diag.flag = FALSE))




