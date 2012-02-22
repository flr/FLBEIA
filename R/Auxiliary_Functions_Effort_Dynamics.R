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
# - capacityRest.SSFB() :: Limits the effort based on total fleet capacity.
# - updateQS.SSFB():: Updates the quota share (QS) in one season and the next 
#          if the quota share does not coincide with the actual catch. 
#          (updates next one only if s < ns and QS(s+1)>0).
#
#
# Dorleta Garcia
# Created: 03/08/2010 10:29:51
# Changed: 01/12/2010 15:23:27  (dga)
# Changed: 2011-02-17 10:31:06  (ssanchez)
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
# capacityRest.SSFB(effs, capacity, order.m)
#       - eff.m[ns,it]  : effort by metier=stock derived from possible catches
#       - capacity[it]  : fleet capacity
#       - order.m[nm,it]: order of importance of the metiers (1st = most important)
#                         (higher importance -> smaller effort reduction)
#                         if not defined, then reduction proportional to effort share
#-------------------------------------------------------------------------------
capacityRest.SSFB <- function(eff.m, capacity, order.m) {
    
    it <- ncol(eff.m)
    mt <- nrow(eff.m)
    
    eff <- apply( eff.m, 2, sum)
    
    eff.share <- eff.m/eff
    
    reduce <- ifelse( eff > capacity, eff - capacity, 0)
    
    for (i in 1:it) {
      if (reduce[i]>0) {
        if ( missing(order.m)) { 
            eff.m[,i] <- eff.m[,i] - reduce[i] * eff.share[,i] 
        } else {
          effs <- eff.m[order.m[,i],i]
          if (length(effs)>=3) {
              red <- Re ( polyroot(c(-capacity[i],sum(effs[-c(1:2)]),effs[2],effs[1]))[1] )
              #red0 <- effRed( eff1=effs[1], eff2=effs[2], effn=sum(effs[-c(1:2)]), capacity=capacity[i])
              eff.m[names(effs[1]),i] <- eff.m[names(effs[1]),i] * red^3
              eff.m[names(effs[2]),i] <- eff.m[names(effs[2]),i] * red^2
              eff.m[names(effs[-c(1:2)]),i] <- eff.m[names(effs[-c(1:2)]),i] * red
              #red; red0; capacity[i]; effs[1]*red^3+effs[2]*red^2+sum(effs[-c(1:2)])*red
          } else if (length(effs)==2) {
              red <- Re ( polyroot(c(-capacity[i],effs[2],effs[1]))[1] )
              eff.m[names(effs[1]),i] <- eff.m[names(effs[1]),i] * red^2
              eff.m[names(effs[2]),i] <- eff.m[names(effs[2]),i] * red
              #red; capacity[i]; effs[1]*red^2+effs[2]*red
          } else {
              red <- Re ( polyroot(c(-capacity[i],effs[1]))[1] )
              eff.m[names(effs[1]),i] <- eff.m[names(effs[1]),i] * red
              #red; capacity[i]; effs[1]*red
          }
        }
      }
    }

  return(eff.m)

}

# effRed: function to estimate the reduction factor to get total.effort=capacity
##################################################################################
#effRed <-  function( eff1, eff2, effn, capacity) {
#
#    totalEff <- function( x, eff1, eff2, effn, capacity) {
#                  return( c( abs( eff1 * x^3 + eff2 * x^2 + effn * x - capacity) ) )
#    }
#
#    res <- optimize( totalEff, c(0,1), tol=.Machine$double.eps^0.25, eff1=eff1, eff2=eff2, effn=effn, capacity)$minimum
#    return(res)
#
#}


#-------------------------------------------------------------------------------
# updateQS.SSFB():: Updates the quota share (QS) in one season and the next 
#          if the quota share does not coincide with the actual catch. 
#          (updates next one only if s < ns and QS(s+1)>0,
#           all the exceeding QS is assigned to the following season).
#       - QS[ns,it]   : quota share (fleet-season) for all the season of certain year.
#       - TAC[it]     : total year TAC for the stock
#       - B[it]       : Seasonal biomass
#       - E[it]       : Predicted effort
#       - efs.m[mt,it]: Effort share of each metier
#       - q.m[mt,it], alpha[mt,it], beta[mt,it] : q and CobDog parameters
#       - season      : the season for which quota share must be updated
#
#   x_s (original quota share)
#   TAC and C (the catch restriction and the actual catch)
#   x_s' = x_s*(C/TAC) (actual quota share)
#           x_i' = x_i + (x_s - x_s'), if x_i>0
#           x_i' = 0                 , if x_i=0
#   x_i' the new quota share for season s+1<=ns
#-------------------------------------------------------------------------------
updateQS.SSFB <- function(QS, TAC, B , E , efs.m, q.m, alpha.m, beta.m, season){

    C    <- colSums(cobbDoug(B = B, E = E, efs.m = efs.m, alpha.m = alpha.m, beta.m = beta.m, q.m = q.m))
    x_s  <- QS[season, ]
    x_s. <- C/TAC
    QS[season,] <- x_s.
    
    if(season == dim(QS)[1]) return (QS)
    
    x_i  <- QS[season+1,]
    x_i. <- ifelse(x_i > 0, x_i  + (x_s - x_s.), x_i)
    x_i. <- ifelse(x_i. < 0, 0, x_i.)
    
    QS[season+1,] <- x_i.
    
    return(QS)
    
}

