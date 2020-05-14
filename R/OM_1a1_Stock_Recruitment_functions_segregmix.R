
#' Stock-Recruitment models in FLBEIA: segrexmix - Mixture of 2 segmented regression stock-recruitment models fit
#' 
#' Model formulation:
#' 
#' \deqn{R = \mathbf{ifelse}(S \leq b, a S, \mathbf{ifelse}(S \leq B, a b, A B uncAdd))}
#'       { R = ifelse(S <= b, a*S, ifelse(S <= B, a*b, A*B*uncAdd))}
#' \emph{a} is the slope of the recruitment for stock levels at or below \emph{b},
#' \eqn{a b}{a*b} is the mean recruitment for stock levels above \emph{b} and at or below \emph{B}, and
#' \eqn{A B uncAdd}{A*B*uncAdd} is the mean recruitment for stock levels above \emph{B}.
#' Where \emph{a, b, A, B} > 0.
#' 
#' Additional stock-recruitment (SR) models to the ones provided in FLCore package.
#' 
#' @author The FLBEIA Team
#' @seealso \link{SRModels}, \linkS4class{FLSR}, \linkS4class{FLModel}


# rec is calculated as:
# - a * SSB       , if     SSB <= b
# - a * b         , if b < SSB <= B
# - A * B * uncAdd, if B < SSB
# 
# where uncAdd corresponds to the uncertainty added when SSB > B,
# with uncAdd = @uncertainty ^ ((sd of SR residuals for SSB <= B) / (sd of SR residuals for SSB <= B) - 1)

# segregmix {{{

#' @name segremix
#' @aliases segregmix

segregmix <- function () 
{
  logl <- function(a, b, A, B, rec, ssb, uncAdd) {
    loglAR1(log(rec), FLQuant(log(ifelse( c(ssb) <= b, a * c(ssb), 
                                          ifelse( c(ssb) <= B, a * b, 
                                                  A * B * uncAdd))), dimnames = dimnames(ssb)))
  }
  
  model <- rec ~ FLQuant(ifelse( c(ssb) <= b, a * c(ssb), 
                                 ifelse( c(ssb) <= B, a * b, 
                                         A * B * uncAdd)), dimnames = dimnames(ssb))
  
  initial <- structure(function(rec, ssb) {
    return(FLPar(a = median(c(rec)/c(ssb), na.rm = TRUE), 
                 b = median(c(ssb), na.rm = TRUE), 
                 A = median(c(rec)/c(ssb), na.rm = TRUE), 
                 B = median(c(ssb), na.rm = TRUE)))
  }, lower = rep(0, 0, 0, 0), upper = rep(Inf, 4))
  
  return(list(logl = logl, model = model, initial = initial))
}

# }}}
