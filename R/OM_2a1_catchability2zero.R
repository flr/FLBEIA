#-------------------------------------------------------------------------------
#       q2zero(fleets, fleets.ctrl, year, season,...)
# 
#   Output: Updated FLFleets 
# 
#   Sets the catchability and or effort share to 0, if possible,
#   when the TAC of a stock is 0:
# 
#       1. Checks if TAC=0 for any stock.
#
#       2. If so updates capturabiblity and effort.share parameters,
#          depending on fleets.ctrl[[fl]]$q2zero values:
#
#       2.1. fleets.ctrl[[fl]]$q2zero[st,mt]=NA (stock st is not caught by metier mt) 
#            --> no changes
#       2.2. fleets.ctrl[[fl]]$q2zero[st,mt]= 0 (stock st is caught by metier mt and cannot be discriminated)
#            --> fleets[[fl]]@metiers[[mt]]@effshare = 0
#       2.3. fleets.ctrl[[fl]]$q2zero[st,mt]= 1 (stock st is caught by metier mt, but can be discriminated) --> 
#            --> fleets[[fl]]@metiers[[mt]]@catches[[st]]@catch.q = 0
#
#
# Sonia Sanchez
# Created: 28/11/2017
# Changed: 01/12/2017 10:52:52
#-------------------------------------------------------------------------------



catchability2zero <- function(fleets, flnm, advice, fleets.ctrl, year) {
  
  it  <- dim(fleets[[1]]@effort)[6]
  ns  <- dim(fleets[[1]]@effort)[4]
  
  mtnms <- names(fleets[[flnm]]@metiers)
  stnms <- catchNames(fleets[[flnm]])
  
  fl.sel <- fleets.ctrl[[flnm]]$q2zero
  
  yr.share <- sapply( stnms, function(x) advice$quota.share[[x]][flnm,year,,,,, drop=T]) # [nst,ni]
  ss.share <- lapply( stnms, function(x) fleets.ctrl$seasonal.share[[x]][flnm,year,,,,, drop=T]) # list[nst] [ns,ni]
  
  if(all(advice$TAC[,year,] > 0) & all(yr.share>0) & all(unlist(ss.share)>0)) return(fleets)
  
  # temp_q <- temp_efs <- numeric(it) + 1 # 1: if unchanged, or 0: if need to convert to 0
  # use sweep function or transform to FLQuants
  
  for(i in 1:it){
    
    TAC <- advice$TAC[,year,,,,i,drop=T]
    
    yr.share <- sapply( stnms, function(x) advice$quota.share[[x]][flnm,year,,,,i, drop=T]) # [nst]
    ss.share <- sapply( stnms, function(x) fleets.ctrl$seasonal.share[[x]][flnm,year,,,,i, drop=T]) # [ns,nst]
    
    if(all(TAC > 0) & all(yr.share>0) & all(ss.share>0))  next
    
    for(st in catchNames(fleets[[flnm]])){

      if (TAC[st]==0 | yr.share[st]==0) {
        
        mtst <- unlist(lapply(mtnms, function(x) if(st %in% names(fleets[[flnm]]@metiers[[x]]@catches)) x))
        
        for (mt in mtst) {
          if (is.na(fl.sel[st,mt])) {
            next()
          } else if (fl.sel[st,mt]==1) { # if TACstk=0 + st can be discriminated in mt --> catch.q[mt_st]=0 
            fleets[[flnm]]@metiers[[mt]]@catches[[st]]@catch.q[,year,,,,i] <- 0
          } else if (fl.sel[st,mt]==0) { # if TACstk=0 + not possible to discriminate st in mt --> effshare[mt]=0
            efs.m <- sapply(mtnms, function(x) fleets[[flnm]]@metiers[[x]]@effshare[,year,,,,i,drop=T])
            for (m in mtst[!(mtst %in% mt)])
              fleets[[flnm]]@metiers[[mt]]@effshare[,year,,,,i] <- efs.m[m] + efs.m[m] * efs.m[mt]/sum(efs.m[names(efs.m)!=mt])
            fleets[[flnm]]@metiers[[mt]]@effshare[,year,,,,i] <- 0
            for (stk in names(fleets[[flnm]]@metiers[[mt]]@catches))
              fleets[[flnm]]@metiers[[mt]]@catches[[stk]]@catch.q[,year,,,,i] <- 0
          }
        }
      
      } else if (ns>1) {
        if (any(ss.share[,st]==0)) {
          for (ss in 1:ns) if (ss.share[ss,st]==0)
            for (mt in mtst) {
              if (is.na(fl.sel[st,mt])) {
                next()
              } else if (fl.sel[st,mt]==1) { # if TACstk=0 + st can be discriminated in mt --> catch.q[mt_st]=0 
                fleets[[flnm]]@metiers[[mt]]@catches[[st]]@catch.q[,year,,ss,,i] <- 0
              } else if (fl.sel[st,mt]==0) { # if TACstk=0 + not possible to discriminate st in mt --> effshare[mt]=0
                efs.m <- sapply(mtnms, function(x) fleets[[flnm]]@metiers[[x]]@effshare[,year,,ss,,i,drop=T])
                for (m in mtst[!(mtst %in% mt)])
                  fleets[[flnm]]@metiers[[mt]]@effshare[,year,,ss,,i] <- efs.m[m] + efs.m[m] * efs.m[mt]/sum(efs.m[names(efs.m)!=mt])
                fleets[[flnm]]@metiers[[mt]]@effshare[,year,,ss,,i] <- 0
                for (stk in names(fleets[[flnm]]@metiers[[mt]]@catches))
                  fleets[[flnm]]@metiers[[mt]]@catches[[stk]]@catch.q[,year,,ss,,i] <- 0
              }
            }
        }
      }

    }
  }

  return(fleets)
  
}

