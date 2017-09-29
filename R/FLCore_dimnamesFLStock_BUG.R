# FLStock.R - FLStock class and methods
# FLCore/R/FLStock.R

# CHANGES added by FLBEIA team to avoid a bug, changes to be incorporated in next FLCore release

# Copyright 2003-2016 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, EC JRC


# dimnames<- {{{
setMethod('dimnames<-', signature(x='FLStock', value='list'),
          function(x, value)
          {
            slots <- getSlotNamesClass(x, 'FLQuant')
            aslots <- c('catch', 'landings', 'discards', 'stock')
            for(i in slots[!slots %in% aslots])
              dimnames(slot(x, i)) <- value
            
            # range
            vnames <- names(value)
            if('year' %in% vnames)
              range(x, c('minyear','maxyear')) <- value[['year']][c(1, length(value[['year']]))]
            if(dims(x)$quant %in% vnames) {
              # BUG What if age='all'
              range(x, c('min','max', 'plusgroup')) <- suppressWarnings(as.numeric(
                value[[dims(x)$quant]][c(1, rep(length(value[[dims(x)$quant]])),2)]))
            }
            
            value <- value[names(value) != dims(x)$quant]
            if(length(value) > 0)
            {
              for (i in aslots)
                dimnames(slot(x, i)) <- value
            }
            
            return(x)
          }
) # }}}
