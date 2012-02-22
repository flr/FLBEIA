#-------------------------------------------------------------------------------
#                       LOW LEVEL FUNCTIONS
#
# Functions that operate with non FLR objects and that perform
# straightforward calculations.
#
# BIOLOGICAL.
# - pope(nay, may, cay, ...):     nay, may, cay ~ matrix[na,it]
# - exponential(nay, may, fay, ...):  nay, may, fay ~ matrix[na,it]
# 
# Dorleta Garcia
# Created: 03/08/2010 10:29:51
# Changed: 01/12/2010 15:22:23  (dga)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# pope() :: project the pop. 1 step forward using Pope's approximation.
#-------------------------------------------------------------------------------
pope <- function(nay, may, cay, ...){
            return((nay*exp(-may/2) -  cay)*exp(-may/2))
            }

#-------------------------------------------------------------------------------
# exponential() :: - project the pop. 1 step forward using Exponential Survival.
#                  - Output: Matrix[na,it]
#-------------------------------------------------------------------------------
exponential <- function(nay, may, fay, ...){
            return((nay*exp(- may - fay)))
            }

## replaceRec {{{
## FLStock <- replaceXSA(FLStock,...)
# Replace the cohorts that are born after (ny - ny.rep) with the geometric mean recruiment
# obtained using the recruitment estimated for (ny.as - ny.rep + 1 - ny.use):(ny.as - ny.rep).

replaceRec <- function(stock, ny.rep, ny.use){

    d <- dim(stock@m)

    it <- d[6]

    ny <- dim(stock@stock.n)[2]
    ny1 <- ny - ny.rep + 1      # First replacement year.
    ny0 <- ny1 - ny.use        # First year to compute the mean recruitment.

    f <- unname(unclass(stock@harvest))
    n <- unname(unclass(stock@stock.n))
    m <- unname(unclass(stock@m))
    cn <- unname(unclass(stock@catch.n))

    # Recruitment Numbers and Fishing Mort.

    for(k in 1:it){
        nrec <- length(n[1, ny0:(ny0 + ny.use - 1),,,,k])
        n[1,ny1:ny,,,,k] <- (prod(n[1, ny0:(ny0 + ny.use - 1),,,,k]))^(1/nrec)


        for(y in ny1:ny){
            f[1,y,,,,k] <- findF(m[1,y,,,,k, drop = F], n[1,y,,,,k, drop = F], cn[1,y,,,,k, drop = F])
        }

        if(ny.rep > 1){
            for(i in 2:ny.rep){
                for(y in (ny1+i-1):ny){
                    n[i,y,,,,k] <- n[i-1,y-1,,,,k]*exp(-(f[i-1,y-1,,,,k, drop = F] + m[i-1,y-1,,,,k, drop = F]))
                    f[i,y,,,,k] <- findF(m[i,y,,,,k, drop = F], n[i,y,,,,k, drop = F], cn[i,y,,,,k, drop = F])
                }
            }
        }
    }

   stock@harvest[] <- f
   stock@stock.n[] <- n

   return(stock)
}

findF <-  function(Ma, Na, Ca){

    d <- dim(Ma)
    dage <- d[1]
    diter <- d[6]

    M <- array(Ma, dim = d)
    N <- array(Na, dim = d)
    Ca <- array(Ca, dim = dim(Ca))  # CAA has no iterations.

    f <- function(Fa,Ma,Na,Ca) return(c(abs(Ca - Fa/(Fa + Ma)*(1-exp(-(Fa+Ma)))*Na)))

    res <- array(dim = d)

    for(i in 1:diter){
        for(a in 1:dage){
                res[a,,,,,i] <- optimize(f,c(0,10), Ma = Ma[a,,,,,i], Na = Na[a,,,,,i], Ca = Ca[a,,,,,i])$minimum

        }
    }
    return(res)
}
