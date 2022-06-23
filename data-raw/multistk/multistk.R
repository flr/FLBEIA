#-------------------------------------------------------------------------------  
# multistk dataset: 
# code to generate data in multiark.RData
# 
# Created: Sonia Sanchez-Maroño -  2022-06-23
# Changed: )
#------------------------------------------------------------------------------- 

# multistk.r - code to generate data in multistk.RData
# FLBEIA/data-raw/multistk/multistk.r

# Copyright: AZTI, 2022
# Author: Sonia Sanchez-Maroño (AZTI) (<flbeia@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


###############################################################################
# TITLE:        Test Case: 1 stock; 3 season; 2 iterations
#      
# Create FLStock with multiple seasons and iterations for calculating seasonal reference points
###############################################################################

library(FLBEIA)

nit <- 3 # Number of iterations

stkn <- "STK"               # Stock name
ages <- 0:6                 # age classes

nyr <- 5

# FLStock with 2 seasons and 3 iterations
stk <- FLStock( name = stkn, 
                stock.wt = FLQuant(NA, dim = c(length(ages),nyr,1,2,1,nit), 
                                   dimnames = list(age=ages,year=1:nyr,iter=1:nit)))

# nit iterations
stk <- propagate( stk, nit)

units(harvest(stk))  <- "f"
units(stock.wt(stk)) <- "kg"

# Mean weight at age by semester (in kg)
stock.wt(stk)[,,,1,] <- c(0.0000, 0.0092, 0.0263, 0.0372, 0.0425, 0.0448, 0.0458)
stock.wt(stk)[,,,2,] <- c(0.0021, 0.0182, 0.0327, 0.0404, 0.0439, 0.0454, 0.0460)

# F range to estimate mean F
range(stk)[c("minfbar","maxfbar")] <- c(min=1,max=3)

# Prop of spawn before M and F
ssb.ss <- rec.ss <- 2
harvest.spwn(stk)[,,,1,] <- m.spwn(stk)[,,,1,] <- 0.25
harvest.spwn(stk)[,,,2,] <- m.spwn(stk)[,,,2,] <- 0

# Selectivity
harvest(stk) <- c(0, 0.85, 1, 1, 1, 1, 1)

# Natural mortality and maturity
m(stk)     <- c(2.9, 0.75, 0.56, 0.50, 0.48, 0.47, 0.47)
m(stk)[1,,,1,] <- 0    # different for age 0, as rec occurs mid-year (not 1st Jan)
mat(stk) <- c(0, 0.5, 1, 1, 1, 1, 1)

# Initial numbers: virgin biomass
stock.n(stk)[ ,1,,1,] <- c(12778112, 8893838, 1949082, 634041, 231392, 88156, 55933)
stock.n(stk)[ ,1,,2,] <- stock.n(stk)[,1,,1,] * exp(-m(stk)[,1,,1,])
stock.n(stk)[1,1,,1,] <- 0 # recruitment occurs in the 2nd semester
units(stock.n(stk)) <- "1000"
stock(stk) <- quantSums(stock.n(stk) * stock.wt(stk))

# Catches: 0 in the 1st simulation year
catch.wt(stk)[,,,1,] <- landings.wt(stk)[,,,1,] <- discards.wt(stk)[,,,1,] <- c(0.0004, 0.0137, 0.0297, 0.0390, 0.0433, 0.0451, 0.0459)
catch.wt(stk)[,,,2,] <- landings.wt(stk)[,,,2,] <- discards.wt(stk)[,,,2,] <- c(0.0052, 0.0225, 0.0352, 0.0415, 0.0444, 0.0456, 0.0461)


#==============================================================================
# SAVE
#==============================================================================

multistk <- stk

save( multistk, file = file.path("data", "multistk.RData"), compress="xz")

