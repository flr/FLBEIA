# Model Free HCRs used for Hake in the Gadoid Simposium.

#-------------------------------------------------------------------------------
# Implementation of the HCR presented in Pomarede 2010. 
#-------------------------------------------------------------------------------
#' @rdname annualTAC
pidHCR <- function(indices, advice, advice.ctrl, year, stknm,...){
  
  Idnm <- advice.ctrl[[stknm]][['index']]  # either the name or the position of the index in FLIndices object.
  Id <- indices[[stknm]][[Idnm]]@index
  
  # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
  year.or <- year
  yrnm    <- dimnames(advice$TAC)[[2]][year]
  year    <- which(yrnm == dimnames(Id)[[2]])
  
  Kp  <- advice.ctrl[[stknm]][['ref.pts']]['Kp',]
  Ki  <- advice.ctrl[[stknm]][['ref.pts']]['Ki',]
  Kd  <- advice.ctrl[[stknm]][['ref.pts']]['Kd',]
  tau <- advice.ctrl[[stknm]][['ref.pts']]['tau',]
  
  alpha  <- advice.ctrl[[stknm]][['ref.pts']]['alpha',]
  
  ey <- numeric(tau)
  
  # Calculate ey, the divergence between the index and the reference level.
  for(y in 0:tau){
    Bnow <- (yearSums(Id[,(year-y-1):(year-y-2)])/2)[drop=T] # [it]
    Bref <- (yearSums(Id[,(year-y-3):(year-y-5)])/3)[drop=T] # [it]
    Brat <- Bnow/Bref  [drop=T] # [it]
    ey[y+1] <- log(Brat)
  }
  
  # Calculate the control signal.
  uy <- exp( Kp*ey[1] + Ki*sum(ey) + Kd*(ey[1]-ey[2]))
  
  if(!is.numeric(alpha)) alpha <- Inf # no constraint in TAC variation.
  
  # Apply +/- alpha% TAC constraint.
  cat("..........", uy, "............\n")
  if( uy > (1 + alpha)) uy <- 1 + alpha
  else if(uy < (1 - alpha)) uy <- 1 -  alpha  
  
  
  advice$TAC[stknm,year.or+1,] <- advice$TAC[stknm,year.or,]*uy
  
  return(advice)
}



#-------------------------------------------------------------------------------
# Implementation of the HCR presented in Pomarede 2010. 
#-------------------------------------------------------------------------------
#' @rdname annualTAC
pidHCRItarg <- function(indices, advice, advice.ctrl, year, stknm,...){
  
  Idnm <- advice.ctrl[[stknm]][['index']]  # either the name or the position of the index in FLIndices object.
  Id <- indices[[stknm]][[Idnm]]@index
  
  # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
  year.or <- year
  yrnm    <- dimnames(advice$TAC)[[2]][year]
  year    <- which(yrnm == dimnames(Id)[[2]])
  
  Kp  <- advice.ctrl[[stknm]][['ref.pts']]['Kp',]
  Ki  <- advice.ctrl[[stknm]][['ref.pts']]['Ki',]
  Kd  <- advice.ctrl[[stknm]][['ref.pts']]['Kd',]
  tau <- advice.ctrl[[stknm]][['ref.pts']]['tau',]
  Itarg <- advice.ctrl[[stknm]][['ref.pts']]['Itarg',]
  
  alpha  <- advice.ctrl[[stknm]][['ref.pts']]['alpha',]
  
  ey <- numeric(tau)
  
  # Calculate ey, the divergence between the index and the reference level.
  for(y in 0:tau){
    Bnow <- (yearSums(Id[,(year-y-1):(year-y-2)])/2)[drop=T] # [it]
    Brat <- (Bnow/Itarg)[drop=T] # [it]
    ey[y+1] <- log(Brat)
  }
  
  # Calculate the control signal.
  uy <- exp( Kp*ey[1] + Ki*sum(ey) + Kd*(ey[1]-ey[2]))
  
  if(!is.numeric(alpha)) alpha <- Inf # no constraint in TAC variation.
  
  # Apply +/- alpha% TAC constraint.
  cat("..........", uy, "............\n")
  if( uy > (1 + alpha)) uy <- 1 + alpha
  else if(uy < (1 - alpha)) uy <- 1 -  alpha  
  
  
  advice$TAC[stknm,year.or+1,] <- advice$TAC[stknm,year.or,]*uy
  
  return(advice)
}


#-------------------------------------------------------------------------------
# Dorleta Garcia
# 09/11/2013 18:44:26
# Implementation of the HCR presented in Little 2011. 
# We have include the contraint min(res,Cmax) not to allow very high catches it can be turn off
# setting Cmax equal to a really high number or Inf.
#-------------------------------------------------------------------------------
#' @rdname annualTAC
little2011HCR <- function(indices, advice, advice.ctrl, year, stknm,...){
  
  Idnm <- advice.ctrl[[stknm]][['index']]  # either the name or the position of the index in FLIndices object.
  Id <- indices[[stknm]][[Idnm]]@index
  
  # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
  year.or <- year
  yrnm    <- dimnames(advice$TAC)[[2]][year]
  year    <- which(yrnm == dimnames(Id)[[2]])
  
  Ctarg  <- advice.ctrl[[stknm]][['ref.pts']]['Ctarg',]
  Ilim   <- advice.ctrl[[stknm]][['ref.pts']]['Ilim',]
  Itarg  <- advice.ctrl[[stknm]][['ref.pts']]['Itarg',]
  Cmax <- advice.ctrl[[stknm]][['ref.pts']]['Cmax',]

  Iy <- (yearSums(Id[,(year-1):(year-2)])/2)[drop=T] # [it]
  
  TAC <- pmin(Ctarg*pmax(0,(Iy-Ilim)/(Itarg-Ilim)),Cmax)
  
  cat("..........", TAC, "............\n")
  
  advice$TAC[stknm,year.or+1,] <- TAC
  
  return(advice)
}

