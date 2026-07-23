#-------------------------------------------------------------------------------  
# TAXATION system related functions
# Created: Sonia Sanchez-Maroño - 2026-07-15
# Changed: 
#------------------------------------------------------------------------------- 

# Auxiliary_Functions_taxation.R - functions to estimate taxes
# ./R/Auxiliary_Functions_taxation.R

# Copyright: AZTI, 2026
# Author: Sonia Sanchez-Maroño (AZTI)
#
# Distributed under the terms of the XXX Public Licence (XXX) V.X.X.


#----------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# convexTax(fleet, fleet.ctrl, advice)
#-------------------------------------------------------------------------------
# Parameters:
# - cat.flst    : fleet's catches of a specific stock.
# - qsh.flst    : fleet's quota share for the stock.
# - tac.st      : stock's TAC.
# - beta        : convex taxation's beta parameter (for this fleet and stock).
# - gamma       : convex taxation's gamma parameter (for this fleet and stock). 
# - tax.rewards : TRUE(default)/FALSE. Indicates if int the convex taxation rewards part should be included or not.

convexTax <- function(cat.flst, qsh.flst, tac.st, 
                      beta = fleet.ctrl[[st]][['beta']], gamma = fleet.ctrl[[st]][['gamma']], 
                      tax.rewards = fleet.ctrl[[st]][['tax.rewards']]) { 
  
  # Given:
  # - Quadratic Tax Instrument: 
  #     T(h,beta,gamma) = beta * h + gamma/2 * h^2
  #   with h=harvest, beta>=0, gamma>=0
  # - For specific fleet and stock (i)
  #   * Taxes    : Ti = si * T(hi/si,beta,gamma)
  #           si = quota share (percent)
  #           hi = individual catch
  #   * Subsidies: R(si) = si * T(TAC,beta,gamma)
  #  If compliance (i.e. hi=TAC*si) --> Taxes - Subsidies = 0 
  #
  # FORMULATION
  # tax.flst = taxes - rewards
  # where:
  #      taxes   = beta * cat.flst + gamma/2 * cat.flst^2/qsh.flst
  #      rewards = beta * tac.st * qsh.flst + gamma/2 * tac.st^2 * qsh.flst
  #
  # used formulation where: Cr.f = QS * tac
  # Taxes should only apply when there is an overshoot, but the undershoot, should be rewarded (i.e. subsidies).
  
  taxes <- beta * cat.flst + gamma/2 * (cat.flst^2 / qsh.flst)
  
  if (tax.rewards == TRUE) {
    rewards <- beta * tac.st * qsh.flst + gamma/2 * (tac.st)^2 * qsh.flst
  } else
    rewards <- 0
  
  tax <- taxes - rewards
  
  return(tax)
  
}

