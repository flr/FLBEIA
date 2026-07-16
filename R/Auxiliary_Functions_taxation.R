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
convexTax <- function(cat.flst, qsh.flst, tac.st, 
                      beta = fleet.ctrl[[st]][['beta']], gamma = fleet.ctrl[[st]][['gamma']]) { 
  
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
  # tax.flst = taxes - rewards = 
  #   = beta * (cat.flst - tac.st * qsh.flst) + 
  #     + gamma/2 * ((cat.flst/qsh.flst)^2 * qsh.flst - tac.flst^2 * qflst)
  # used formulation where: Cr.f = QS * tac
  # Taxes should only apply when there is an overshoot, but the undershoot, should be rewarded (i.e. subsidies).
  
  tax <- beta * (cat.flst - tac.st * qsh.flst) + 
    gamma/2 * (cat.flst^2 / qsh.flst - (tac.st)^2 * qsh.flst)
  
  return(tax)
  
}

