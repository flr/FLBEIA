#-------------------------------------------------------------------------------
#                           Assessments
#   - assessment.mp:
#         update, stock.n, stock and harvest for each stock, 
#   - NoAssessment.
#
# Dorleta Garcia
# Created: 21/12/2010 07:55:33
# Changed: 21/12/2010 07:55:38
# Changes: 2012-06-15 12:59:05  Sonia Sanchez - for allowing assessment in different seasons and multiannual advice
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# assessment.mp(stocks, fleets.obs, indices, assess.ctrl)   
#  For the time being only single stock assessments.
#-------------------------------------------------------------------------------
assessment.mp <- function(stocks, fleets.obs, indices, covars=covars, assess.ctrl, datayr, season, stknm){

    st <- stknm
     
    if(assess.ctrl[[st]][['assess.model']] != "NoAssessment") {
      
        if (assess.ctrl[[st]]$ass.curryr & season==dim(stocks[[1]]@stock.n)[4]) datayr <- datayr+1
        
        # trim the indices, fron index specific initial year to the assessment year.
        indST <- indices[[st]]
        if(!is.null(indST))
            indST <- FLIndices(lapply(indices[[st]],function(x) trim(x, year = dimnames(x@index)[[2]][1]:datayr)))
        
        if(assess.ctrl[[st]]$work_w_Iter == TRUE){  # The assessment model works with iters => all the iterations are 'adjusted' in one run.
            res <- eval(call(assess.ctrl[[st]][['assess.model']], stock = stocks[[st]], indices = indST, control = assess.ctrl[[st]]$control,covars=covars))
            stock.n(stocks[[st]]) <- res$stock@stock.n
            harvest(stocks[[st]]) <- res$stock@harvest
            covars <- res$covars
           # stock(stocks[[st]])   <- res@stock
        }
        else{  # The assessment model _DOES NOT_ work with iters =>  the iterations are 'adjusted' one by one.
            it <- dim(stocks[[st]]@m)[6]
            
            for(i in 1:it){ 
                res <- eval(call(assess.ctrl[[st]][['assess.model']], stock = iter(stocks[[st]],i), 
                            indices = iter(indST,i), control = assess.ctrl[[st]]$control))
                iter(stock.n(stocks[[st]]),i) <- res$stock@stock.n
                iter(harvest(stocks[[st]]),i) <- res$stock@harvest
                covars <- res$covars
            }
        }
        stock(stocks[[st]])   <- computeStock(stocks[[st]])
        units(harvest(stocks[[st]])) <- assess.ctrl[[st]]$harvest.units
    }
    
    return(list(stocks=stocks,covars=covars))

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
FLXSA2flbeia <- function(stock, indices, control,covars=covars) return(list(stock=FLXSA::FLXSA(stock, indices, control, diag.flag = FALSE),covars=covars))




