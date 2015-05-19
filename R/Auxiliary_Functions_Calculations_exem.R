# Dorleta Garcia 2015/01/14
# Functions
#
# overfishing: Are the stock suffering of overfishing, returns a matrix [nst,nit]

#-------------------------------------------------------------------------------------------
# overfishing: Are the stock suffering of overfishing, returns a matrix [nst,nit]
#-------------------------------------------------------------------------------------------
overfishing <- function(biols, fleets, advice.ctrl, year){
  if(is.numeric(year)) year <- dimnames(biols[[1]]@n)[[2]][year]
  obj <- list(biols = biols, fleets = fleets)
  f    <- F_flbeia(obj, years = year)
  
  Fref <- dimnames(x$ref.pts)[[1]]
  Fref <- Fref[which(substr(Fref, 1,1) == 'F')]
  
  fmsy <- lapply(advice.ctrl, function(x) ifelse(is.null(x$ref.pts[Fref,]), NA, x$ref.pts[Fref,]))
  res <- sapply(names(fmsy), function(x) fmsy[[x]] < f[x,,])
  res <- ifelse(is.na(res), FALSE, res)
  if(!is.null(dim(res))) res <- t(res) # if there are iterations put them in columns not in rows.
  
  res <- as.matrix(res, nrow = (length(biols)))
  
  return(res)  
}