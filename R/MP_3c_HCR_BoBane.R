#-------------------------------------------------------------------------------
#   HCRs for the Bay of Biscay anchovy
#   - annualTAC.
#
# Sonia Sánchez
# Created: 04/06/2012 15:55:35
# Changed: 16/10/2013 12:01:18
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# General HCR for BoBane long term management plan - SSB based
# 
#  TAC advice depending on SSB in relation to Btrigger points is:
#    (July-June)
#      - 0               , ssb <= Btrig1
#      - Tacmin          , Btrig1 < ssb <= Btrig2
#      - alpha+gamma*ssb , Btrig2 < ssb <= Btrig3
#      - TACmax         , ssb > Btrig3
#      If not defined alpha=0, Btrig3 = TACmax/gamma
#    (January-December)
#      - 0               ,      0 < (SSB|TACmin) <= Btrig1
#      - Tacmin          , Btrig1 < (SSB|TACmin) <= Btrig2
#      - gamma1          , (SSB|gamma*SSB) <= Btrig2 < (SSB|TACmin)
#      - alpha+gamma*ssb , Btrig2 < (SSB|gamma*SSB) <= Btrig3
#      - gamma2          , (SSB|TACmax) <= Btrig3 < (SSB|gamma*SSB)
#      - TACmax          , Btrig3 < (SSB|TAcmax)
#   Constraints: If not defined alpha=0, Btrig3 = TACmax/gamma
#                Btrig1 <= Btrig2 <= Btrig3
#                TACmin <= TACmax
#                (TACmin-alpha)/Btrig2 <= gamma <= (TACmax-alpha)/Btrig3
#
#  - Proposed rule for BoBane long term management plan (July-June):
#    (Units = thousand tons)
#           - 0       , ssb <= 24
#           - 7       , 24 < ssb <= 33
#           - 0.3*ssb , 33 < ssb <= 110
#           - 33      , ssb > 110
#
#  - Proposed rule by the RAC in 2013 (July-June):
#    (Units = thousand tons)
#           - 0         , ssb <= 24
#           - 7         , 24 < ssb < 33
#           - gamma*ssb , 33 < ssb <=58
#           - 25        , ssb > 58
#
#-------------------------------------------------------------------------------


aneHCRs <- function(indices, advice, advice.ctrl, year, season, stknm,...){

  # Rule designed for cases with 2 seasons (i.e. semesters)
  ns <- 2
  
  Idnm <- advice.ctrl[[stknm]][['index']]  # either the name or the position of the index in FLIndices object.
  Id <- indices[[stknm]][[Idnm]]@index
  
  # Year  => Character, because the year dimension in indices does not coincide with year dimension in biol.
  year.or <- year
  yrnm    <- dimnames(advice$TAC)[[2]][year]
  year    <- which(yrnm == dimnames(Id)[[2]])
  
  iter    <- dim(Id)[6]
  
  ref.pts <- advice.ctrl[[stknm]]$ref.pts
  alpha  <- ref.pts['alpha',] # Same units as biomass
  gamma  <- ref.pts['gamma',] 
  TACmin <- ref.pts['TACmin',]
  TACmax <- ref.pts['TACmax',]
  Btrig1 <- ref.pts['Btrig1',]
  Btrig2 <- ref.pts['Btrig2',]
  Btrig3 <- ref.pts['Btrig3',]
  
  # Complete missing values: 
  alpha  <- ifelse( is.na(alpha), 0, alpha)              # If missing alpha => alpha = 0
  Btrig3 <- ifelse( is.na(Btrig3), (TACmax-alpha)/gamma, Btrig3) # If missing Btrig3 => Btrig3 = TACmax/gamma
  
  # Checking reference points
  if ( sum(Btrig1>Btrig2 | Btrig2>Btrig3)>0 )
    stop("Check reference points for '", stknm, "'. In advice.ctrl[['", stknm, "']]$ref.pts, necessarily Btrig1 <= Btrig2 <= Btrig3")
  if ( sum(TACmin>TACmax)>0 )
    stop("Check reference points for '", stknm, "'. In advice.ctrl[['", stknm, "']]$ref.pts, necessarily TACmin <= TACmax")
  if ( sum(round((TACmin-alpha)/Btrig2,2) > round(gamma,2) | round((TACmax-alpha)/Btrig3,2) < round(gamma,2))>0 )
    stop("Check reference points for '", stknm, "'. In advice.ctrl[['", stknm, "']]$ref.pts, necessarily (TACmin-alpha)/Btrig2 <= gamma <= (TACmax-alpha)/Btrig3") 
      
  if (season == ns) {
    
    # select the appropriate abundance index
    b.idx <- Id[,year,]                             # [2,it]
    B     <- quantSums(b.idx)[drop=T]               # [it]
    BP    <- b.idx[1,drop=T]/B                      # [it]
    # and other parameters
    TACs1.perc <- advice.ctrl[[stknm]]$TACs1.perc   # [1]
    tsurv      <- advice.ctrl[[stknm]]$tsurv        # [1]
    param      <- advice.ctrl[[stknm]]$cbbm.params
    G1         <- param[['G']][1,year,drop=T]       # [it]
    G2         <- param[['G']][2,year,drop=T]       # [it]
    M1         <- param[['M']][1,year,drop=T]       # [it]
    M2         <- param[['M']][2,year,drop=T]       # [it]
    S1         <- param[['S']][1,year,drop=T]       # [it]
    S2         <- param[['S']][2,year,drop=T]       # [it]
        
    # Calculate the TAC.
    if (alpha==0 & gamma==0){
      
      TAC <- 0
      
    } else {
      
      # To estimate TAC short term forecast and optimization needed
      fmin <- stf.Ftac( TAC=rep(TACmin,iter), TACs1.perc=TACs1.perc, B=B, BP=BP, G1=G1, G2=G2, M1=M1, M2=M2, S1=S1, S2=S2)
      fmed <- stf.Fgamma( alpha=alpha, gamma=gamma, TACs1.perc=TACs1.perc, tsurv=tsurv, B=B, BP=BP, G1=G1, G2=G2, M1=M1, M2=M2, S1=S1, S2=S2)
      fmax <- stf.Ftac( TAC=rep(TACmax,iter), TACs1.perc=TACs1.perc, B=B, BP=BP, G1=G1, G2=G2, M1=M1, M2=M2, S1=S1, S2=S2)
      
      ssb.min <- ssb.cbbm( B=B, BP=BP, G1=G1, G2=G2, M1=M1, M2=M2, f=fmin, S1=S1, S2=S2, tsurv=tsurv)
      ssb.med <- ssb.cbbm( B=B, BP=BP, G1=G1, G2=G2, M1=M1, M2=M2, f=fmed, S1=S1, S2=S2, tsurv=tsurv)
      ssb.max <- ssb.cbbm( B=B, BP=BP, G1=G1, G2=G2, M1=M1, M2=M2, f=fmax, S1=S1, S2=S2, tsurv=tsurv)
      
      #  Calculate where we are in relation to reference biomasses
      Brefs <- c(0,Btrig1,Btrig2,Btrig3)
      
      # Find where the SSB is in relation to reference points.
      bmin.pos <- apply( matrix(1:iter,1,iter), 2, function(i) findInt(ssb.min[i], Brefs))  # [it]
      bmed.pos <- apply( matrix(1:iter,1,iter), 2, function(i) findInt(ssb.med[i], Brefs))  # [it]
      bmax.pos <- apply( matrix(1:iter,1,iter), 2, function(i) findInt(ssb.max[i], Brefs))  # [it]
      
      for (i in 1:iter)
        if ( bmin.pos[i]==2 & bmed.pos[i]==3 & bmax.pos[i]==4 ) stop("Not able to find a TAC for '", stknm, "' stock")

      
      TAC   <- ifelse(bmin.pos == 1, 0,                                               #      0 < (SSB|TACmin) <= Btrig1
                      ifelse(bmin.pos == 2, TACmin,                                   # Btrig1 < (SSB|TACmin) <= Btrig2
                             ifelse(bmin.pos >= 3 & bmed.pos <= 2, 
                                      stf.TAC(SSB.obj=Btrig2,TACs1.perc=TACs1.perc,B=B,BP=BP,G1=G1,G2=G2,M1=M1,M2=M2,S1=S1,S2=S2,tsurv=tsurv), 
                                                                                        # (SSB|gamma*SSB) <= Btrig2 < (SSB|TACmin)
                                    ifelse(bmed.pos == 3, alpha+gamma*ssb.med,              # Btrig2 < (SSB|gamma*SSB) <= Btrig3
                                           ifelse(bmax.pos == 3 & bmed.pos == 4, 
                                                  stf.TAC(SSB.obj=Btrig3,TACs1.perc=TACs1.perc,B=B,BP=BP,G1=G1,G2=G2,M1=M1,M2=M2,S1=S1,S2=S2,tsurv=tsurv), 
                                                                                        # (SSB|TACmax) <= Btrig3 < (SSB|gamma*SSB)
                                                  TACmax)))))                         # Btrig3 < (SSB|TAcmax)
    }
    
    advice[['TAC']][stknm,year.or+1,,,,] <- TAC
    

  } else {
    
    # select the appropriate abundance index
    b.idx <- Id[,year,drop=T] # [it]
    
    #  Calcuate where we are in relation to reference biomasses
    Brefs <- c(0,Btrig1,Btrig2,Btrig3)
    
    # Find where the SSB is in relation to reference points.
    b.pos <- apply( matrix(1:iter,1,iter), 2, function(i) findInt(b.idx[i], Brefs))  # [it]
    
    # Calculate the TAC.
    if (alpha==0 & gamma==0){
      TAC <- 0
    } else {
      TAC   <- ifelse(b.pos == 1, 0,                         #      0 < (SSB|TACmin) <= Btrig1
                      ifelse(b.pos == 2, TACmin,             # Btrig1 < (SSB|TACmin) <= Btrig2
                             ifelse(b.pos == 3, alpha+gamma*b.idx, # Btrig2 < (SSB|gamma*SSB) <= Btrig3
                                    TACmax)))                # Btrig3 < (SSB|TAcmax)
    }
    
    advice[['TAC']][stknm,year.or,,,,] <- TAC
    
  }
  
  return(advice)

}

#-------------------------------------------------------------------------------
# stf.Fgamma(gamma,TACs1.perc,BP,G1,G2,M1,M2,S1,S2,tsurv) :: estimates the f to meet the condition TAC=gamma*SSB
#-------------------------------------------------------------------------------
# Function to calculate TAC (TAC=alpha+gamma*SSB), given:
# alpha     : intecept from the HCR                                                  [1]
# gamma     : harvest rate from the HCR                                              [1]
# TACs1.perc: percentage of the TAC which is assumed to be taken in the 1st semester [1]
# B         : biomass 1+ (1st January)                                               [it]
# BP        : percentage of age 1                                                    [it]
# Gi        : growth rate at i, i=1,2                                                [it]
# Mi        : mortality rate at age i                                                [it]
# Si        : selectivity at age i, in the 1st season                                [it]
# tsurv     : time of the survey                                                     [1]
stf.Fgamma <- function( alpha, gamma, TACs1.perc, B, BP, G1, G2, M1, M2, S1, S2, tsurv) {
  
  nit <- length(BP) 
  
  of <- function(f,alpha,gamma,TACs1.perc,B,BP,G1,G2,M1,M2,S1,S2,tsurv) 
    return(c(abs( alpha * TACs1.perc / B +
                  BP * (TACs1.perc*gamma*exp((G1-M1-f*S1)*tsurv)-(1-exp((G1-M1-f*S1)*0.5))*((f*S1)/(-G1+M1+f*S1))) +
               (1-BP) * (TACs1.perc*gamma*exp((G2-M2-f*S2)*tsurv)-(1-exp((G2-M2-f*S2)*0.5))*((f*S2)/(-G2+M2+f*S2))))))

  res <- rep(NA,nit)
  
  for (i in 1:nit) {
    res.opt <- optimize(of, c(0,10), alpha=alpha, gamma=gamma, TACs1.perc=TACs1.perc, B=B[i], BP=BP[i], G1=G1[i], G2=G2[i], M1=M1[i], M2=M2[i], S1=S1[i], S2=S2[i], tsurv=tsurv,
                        tol = .Machine$double.eps^0.5)
    if ( res.opt$obj > 0.001 ) warning("The optimizer in the function 'stf.Fgamma' may not have reach a minimum. OBJ=",res.opt$obj )
    res[i]  <- res.opt$min
  }
  
  return(res)
    
}

#-------------------------------------------------------------------------------
# stf.Ftac(TAC,TACs1.perc,B,BP,G1,G2,M1,M2,S1,S2) :: estimates the f given a TAC
#-------------------------------------------------------------------------------
# Function to calculate TAC (TAC=gamma*SSB), given:
# TAC       : TAC value                                                              [1]
# TACs1.perc: percentage of the TAC which is assumed to be taken in the 1st semester [1]
# B         : biomass 1+ (1st January)                                               [it]
# BP        : percentage of age 1                                                    [it]
# Gi        : growth rate at i, i=1,2                                                [it]
# Mi        : mortality rate at age i                                                [it]
# Si        : selectivity at age i, in the 1st season                                [it]
stf.Ftac <- function( TAC, TACs1.perc, B, BP, G1, G2, M1, M2, S1, S2) {
  
  nit <- length(BP) 
  
  of <- function(f,TAC,TACs1.perc,B,BP,G1,G2,M1,M2,S1,S2) 
    return(c(abs( TAC * TACs1.perc - B * BP * (1-exp((G1-M1-f*S1)*0.5)) * ((f*S1)/(-G1+M1+f*S1)) - 
                                  B * (1-BP) * (1-exp((G2-M2-f*S2)*0.5)) * ((f*S2)/(-G2+M2+f*S2)) )))
  
  res <- rep(NA,nit)
  
  for (i in 1:nit) {
    res.opt <- optimize(of, c(0,10), TAC=TAC[i], TACs1.perc=TACs1.perc, B=B[i], BP=BP[i], G1=G1[i], G2=G2[i], M1=M1[i], M2=M2[i], S1=S1[i], S2=S2[i],
                        tol = .Machine$double.eps^0.5)
    if ( res.opt$obj > 0.001 ) warning("The optimizer in the function 'stf.Ftac' may not have reach a minimum. OBJ=",res.opt$obj )
    res[i]  <- res.opt$min
  }
  
  return(res)
  
}

#-------------------------------------------------------------------------------
# stf.TAC(SSB.obj,TAC,TACs1.perc,B,BP,G1,G2,M1,M2,S1,S2, tsurv) :: estimates the TAC given a desired SSB
#-------------------------------------------------------------------------------
# Function to calculate TAC, given:
# SSB.obj   : desired SSB                                                            [1]
# TACs1.perc: percentage of the TAC which is assumed to be taken in the 1st semester [1]
# B         : biomass 1+ (1st January)                                               [it]
# BP        : percentage of age 1                                                    [it]
# Gi        : growth rate at i, i=1,2                                                [it]
# Mi        : mortality rate at age i                                                [it]
# Si        : selectivity at age i, in the 1st season                                [it]
stf.TAC <- function( SSB.obj, TACs1.perc, B, BP, G1, G2, M1, M2, S1, S2, tsurv) {
  
  nit <- length(BP) 
  
  SSB <- SSB.obj
  
  of <- function(f,SSB,B,BP,G1,G2,M1,M2,S1,S2,tsurv) 
    return(c(abs( SSB - B * BP * exp((G1-M1-f*S1)*tsurv) - B * (1-BP) * exp((G2-M2-f*S2)*tsurv) )))
  
  fval <- c()
  
  for (i in 1:nit) {
    fval.opt <- optimize(of, c(0,10), SSB=SSB[i], B=B[i], BP=BP[i], G1=G1[i], G2=G2[i], M1=M1[i], M2=M2[i], S1=S1[i], S2=S2[i], tsurv=tsurv,
                         tol = .Machine$double.eps^0.5)
    if ( fval.opt$obj > 0.001 ) warning("The optimizer in the function 'stf.gamma' may not have reach a minimum. OBJ=",res.opt$obj )
    fval[i]  <- fval.opt$min
  }
  
  C <- ( B * BP * (1-exp((G1-M1-fval*S1)*0.5)) * ((fval*S1)/(-G1+M1+fval*S1)) + 
         B * (1-BP) * (1-exp((G2-M2-fval*S2)*0.5)) * ((fval*S2)/(-G2+M2+fval*S2)) )

  TAC <- C/TACs1.perc

  return(TAC)
  
}

# Function to calculate SSB, when fishing mortality is known:
# B    : biomass 1+ (1st January)
# BP   : percentage of age 1 (in mass)
# Ga   : growth rate at age
# Ma   : mortality rate at age
# f    : fishing mortality rate
# Sa   : selectivity in the 1st season at age
# tsurv: time of the survey
ssb.cbbm <- function( B, BP, G1, G2, M1, M2, f, S1, S2, tsurv) {
  ssb <- B * ( BP * exp((G1-M1-f*S1)*tsurv) + (1-BP) * exp((G2-M2-f*S2)*tsurv) )
  return(ssb)
}

# Find the indices of 'x' in 'vec' (similar to findInterval)
#  - findInt     : v[i(j)]< x[j]<=v[i(j)+1]
#  - findInterval: v[i(j)]<=x[j]< v[i(j)+1]
findInt <- function(x,vec) return(apply(outer(x, vec, ">"),1,sum))

