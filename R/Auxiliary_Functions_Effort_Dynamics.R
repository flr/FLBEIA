#-------------------------------------------------------------------------------
#                       LOW LEVEL FUNCTIONS
#
# Functions that operate with non FLR objects and that perform
# straightforward calculations.
#
# FLEETS.
# - cobbDoug()
# - cobbDougInv()
# - NomEff.CatchRest() :: Calculates the effort per fleet under Quota share = catch restriction
# 
# ** SMFB functions. The following functions have been developed to be used within 
#           SMFB function, but could be useful for other functions.
# - effRest.SMFB() :: Calculates next year effort applying the 'rule' in effort
#                     restriction. It is designed to work within SMFB function.
# - capacityRest.SMFB() :: Limit the effort based on capacity.
# - updateQS.SMFB():: Updates the quota share (QS) in one season and the next 
#          if the quota share does not coincide with the actual catch. 
#          (updates next one only if s < ns).
#
# ** SSFB functions. The following functions have been developed to be used within 
#           SSFB function, but could be useful for other functions.
# - effrule.SSFB() :: Limits the effort based on total fleet capacity and distributes the 
#                     effort between metiers considering expected effort and available quotas.
# - remQ.SSFB():: Calculates fleets' remanent quota of an stock.
#
#
# Dorleta Garcia
# Created: 03/08/2010 10:29:51
# Changed: 01/12/2010 15:23:27  (dga)
# Changed: 2011-02-17 10:31:06  (ssanchez) - new functions added
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# cobbDoug() :: Cobb-Douglass production fun. at single fleet single stock level.
#               Input: * B,E vector[it]
#                      * efs, alphas, betas, qs: matrix[nmetier,it]
#               Output: Matrix[nmetier,it]  (catch per metier)
#-------------------------------------------------------------------------------
cobbDoug <- function(B, E, efs.m, alpha.m, beta.m, q.m, ...){

            Dim <- dim(efs.m)
            E.m <- matrix(E, Dim[1], Dim[2], byrow = T)*efs.m
            B.m <- matrix(B, Dim[1], Dim[2], byrow = T) 

            return(q.m*E.m^alpha.m*B^beta.m)
}


#-------------------------------------------------------------------------------
# cobbDougInv() :: Cobb-Douglass production fun. at single fleet single stock level.
#                  to estimate necessary effort to allow certain level of catches
#                  restriction: only one metier per stock in each fleet
#               Input: * Cr, B, q.m, alpha.m, beta.m: numeric
#               Output: numeric = effort (needed to allow Cr catch level)
#-------------------------------------------------------------------------------
cobbDougInv <- function(Cr, B,  q.m, alpha.m, beta.m,...){
             (Cr *B^-beta.m / q.m) ^ (1/alpha.m)
            }




#-------------------------------------------------------------------------------
# effRule.SMFB(effs, prev.eff, rule)
#       - effs[ns,it]:  Forecasted effort.
#       - prev.eff[it]: Previous year effort.
#       - rule[character]: rule to be applied to select effort from effs.
#-------------------------------------------------------------------------------

effRule.SMFB <- function(effs, prev.eff, rule){
    sts <- rownames(effs)

    if(dim(effs)[1] == 1)                  return(effs)  # The fleet only catch one stock.
    if(rule %in% sts)                     return(effs[rule,])
    if(rule %in% c('max', 'min', 'mean')) return(apply(effs, 2, rule))
    #=> rule = previous.
   # cat('~~~~~~~~~~~~~~ Rule: Previous Effort ~~~~~~~~~~\n')
    eff <- numeric(dim(effs)[2])
    inds <- apply(abs(1-sweep(effs, 2, prev.eff[drop=T], "/")),2, which.min)
    for(i in 1:dim(effs)[2])  eff[i] <- effs[inds[i],i]
    return(eff)
}


#-------------------------------------------------------------------------------
# capacityRest.SMFB(eff, capacity)
#       - eff[it]     :    Output of  effRule.SMFB.
#       - prev.eff[it]:    Previous year effort.
#       - rule[character]: rule to be applied to select effort from effs.
#-------------------------------------------------------------------------------

capacityRest.SMFB <- function(eff, capacity){
    eff <- ifelse(eff > capacity, capacity, eff)
  return(eff)
}


#-------------------------------------------------------------------------------
# - updateQS.SMFB():: Updates the quota share (QS) in one season and the next 
#          if the quota share does not coincide with the actual catch. 
#          (updates next one only if s < ns).
#   * QS[ns,it]: quota share (fleet-season) for all the season of certain year.
#   * Cr[it]: The catch restriction for the season
#   * B[it]: Seasonal biomass.
#   * E[it]: Predicted effort.
#   * efs.m[mt,it], q.m[mt,it], alpha[mt,it], beta[mt,it]: Effort share, q and CobDog parameters.
#   * season: the season for which quota hare mus be updated.
#
#   x_s (original quota share)
#   Cr and C (the catch restriction and the actual catch)
#   x_s' = x_s*(C/Cr) (actual quota share)
#           x_i' = x_i + (x_s - x_s')/(numb. of seas : i>s)
#   x_i' the new quota share for season i such that i>s>=ns.
#-------------------------------------------------------------------------------
updateQS.SMFB <- function(QS, TAC, catch, season){

    x_s  <- QS[season, ]
    x_s. <- catch/TAC
    QS[season,] <- x_s.
    
    if(season == dim(QS)[1]) return (QS)
    
    x_i  <- matrix(QS[-(1:season),], dim(QS)[1]-season, dim(QS)[2])
    x_i. <- x_i  + matrix((x_s - x_s.)*x_i/sum(x_i), dim(x_i)[1], dim(x_i)[2], byrow = TRUE)
    x_i. <- ifelse(x_i. < 0, 0, x_i.)
    
    QS[-(1:season),] <- x_i.
    
    return(QS)
    
}


#-------------------------------------------------------------------------------
# effRule.SSFB(effs, prev.eff, rule)
#       - ef.m[nmt,it]: effort by metier
#       - efs.m[nmt,it]: effort share by metier
#       - rule[character]: rule to be applied to reasign effort if TAC exhausted
#-------------------------------------------------------------------------------
# Function to estimate fleet effort depending on pre-defined effort share and
# available TAC
# Output: ef.m[nmt,it] - forecasted effort by metier

effRule.SSFB <- function( Ba, B, ef.m, efs.m, q.m, alpha.m, beta.m, fleet, fleet.ctrl, flinfo, Cr.f, Er.f){
  
  it      <- ncol(ef.m)
  sts     <- catchNames(fleet)
  
  ef.fl <- apply(ef.m,2,sum)
  
  stock.ctrl <- matrix(1,length(sts), it, dimnames = list(sts, 1:it)) # 0 if stock's TAC exhausted [nst,it]
  for (st in sts) stock.ctrl[st,] <- ifelse(Cr.f[st,] == 0 & ef.m[flinfo[[st]][2],]==0, 0, stock.ctrl[st,])
    
  ef.tr <- rep( 0, it) # fleet effort to be reasigned
  
  for ( i in 1:it) if (ef.fl[i]!=0) {
    
    while ( sum(stock.ctrl[,i]==1)>0 ) {
      
      for(st in sts){
        
        if (stock.ctrl[st,i]<=0) next
        
        mtst <- flinfo[[st]][2]
        
        stock.ctrl[st,i] <- -1
        
        # Expected catches
        catchFun <- paste(fleet.ctrl[[st]][['catch.model']], 'CatchFleet', sep = ".")
        catch <- eval(call(catchFun, Ba = Ba[[st]], B = B[st,], effort = ef.fl, efs.m = efs.m, q.m = q.m[[st]], alpha.m = alpha.m[[st]], beta.m = beta.m[[st]]))
        #catch <- colSums( cobbDoug(Ba = Ba[[st]], B = B[st,], E = ef.fl, efs.m = efs.m, q.m = q.m[[st]][,,,,drop=T], alpha.m = alpha.m[[st]][,,,,drop=T], beta.m = beta.m[[st]][,,,,drop=T]))
        
        # Catch restrictions (i.e. remaining quota)
        if ( catch[i] > Cr.f[st,i] ) {
          
          catch[i] <- Cr.f[st,i]
          eff.new  <- Er.f[mtst,i]
          ef.tr[i] <- ef.tr[i] + (ef.m[mtst,i] - eff.new)
          ef.m[mtst,i]  <- eff.new
          efs.m[mtst,i] <- ef.m[mtst,i]/ef.fl[i]
          stock.ctrl[st,i] <- 0
          
          # Re-allocation of effort
          rule <- fleet.ctrl$effort.realloc
          st.a <- sts[!(sts %in% st)]
          st.a <- st.a[stock.ctrl[st.a,i]!=0]
          fl.a <- sapply( st.a, function(x) flinfo[[x]][2])
          
          if (length(fl.a)==0) next
          
          stock.ctrl[st.a,i] <- 1
          
          if ( is.null(rule) ) { # Reallocate same proportion for all metiers
            ef.m[fl.a,i] <- ef.m[fl.a,i] + ef.tr[i]/length(fl.a)
            efs.m[fl.a,i] <- ef.m[fl.a,i]/ef.fl[i]
            ef.tr[i] <- 0
          } else if ( rule == 'curr.eff') { # Reallocate proportionally to currently allocated effort
            ef.m[fl.a,i] <- ef.m[fl.a,i] + ef.tr[i]*ef.m[fl.a,i]/sum(ef.m[fl.a,i])
            efs.m[fl.a,i] <- ef.m[fl.a,i]/ef.fl[i]
            ef.tr[i] <- 0
          } else { stop( paste("Effort reallocation rule: '", rule,"' is not implemented in SSFB",sep=""))}
        }
        
      } # END sts LOOP
      
    } # END while
    
    #! CHECKING:
    if (ef.tr[i]==0) {
      if ( round(ef.fl[i] - sum(ef.m[, i]),10) != 0) stop('Problem in effRule.SSFB function - efforts not correctly reallocated')
      if ( round(sum( efs.m[,i] - ef.m[,i]/ef.fl[i]),10) != 0) stop('Problem in effRule.SSFB function - effort share not correctly re-estimated')
    }#! end CHECKING
    
  } # END it LOOP
  
  return(ef.m)
}


#-------------------------------------------------------------------------------
# remQ.SSFB (fleets, TAC, QS, ass.ss, year, season, flnm, stknm, ns)
#       - TAC[st,it]:  TAC set in year
#       - QS[[st]][fl,it]: quota share for each stock by fleet
#-------------------------------------------------------------------------------
# Function to estimate the fleet's remanent quota for an stock
# Output: remQ[it]
remQ.SSFB <- function( fleets, TAC, QS, ass.ss, year, season, flnm, stknm, ns){
  
  fleet <- FLFleetsExt(fleets[[flnm]])
  
  it    <- ncol(TAC)
  
  if (is.null(ass.ss)) { ass.ss <- ns } else if (is.na(ass.ss)) { ass.ss <- ns }
  if (!(ass.ss %in% (1:ns))) stop("Assessment season for: '", stknm, "' outside season range in the objects")
  
  TAC.fl <- TAC[stknm,] * QS[[stknm]][flnm,]
  pc <- rep(0,it)
  
  if (stknm %in% catchNames(fleet)) {
    if (ass.ss==ns & season!=1) { # stocks assessed at the end of the year
      pc <- apply(catchWStock( fleet,stknm)[,year,,1:(season-1),],c(2,6),sum)[drop=T]
    } else if (ass.ss!=ns) {
      if (season==1) {
        pc <- apply(catchWStock( fleet,stknm)[,year-1,,(ass.ss+1):ns,],c(2,6),sum)[drop=T]
      } else if (season<=ass.ss) {
        pc <- apply(catchWStock( fleet,stknm)[,year-1,,(ass.ss+1):ns,],c(2,6),sum)[drop=T]+
          apply(catchWStock( fleet,stknm)[,year,,1:(season-1),],c(2,6),sum)[drop=T]
      } else if (season>ass.ss+1) {
        pc <- apply(catchWStock( fleet,stknm)[,year,,(ass.ss+1):(season-1),],c(2,6),sum)[drop=T]
      }
    }    
  }
  
  remQ <- TAC.fl - pc
  remQ <- ifelse( remQ<0, 0, remQ)
  
  return(remQ)
  
} # END function: remQ.SSFB
