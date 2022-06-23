################################################################################
#  Reference points estimation with seasonal steps                             #
#                                                                              #
#------------------------------------------------------------------------------#
#   Sonia Sanchez & Leire Ibaibarriaga (AZTI-Tecnalia)                         #
#   created:  28/11/2019                                                       #  
################################################################################

# Copyright: AZTI, 2019
# Authors: Sonia Sanchez (AZTI)
#          Leire Ibaibarriaga (AZTI)
#
# Distributed under the terms of the GNU GPLv3

#' @title Reference points for a seasonal model
#' 
#' @description Function to estimate reference points when having a seasonal model
#'
#' @param stk An FLStock object.
#' @param B0 The value of the virgin biomass. 
#' @param R0 The expected recrutiment in the virgin population.
#' @param rec.ss The recruitment season (numeric). Default value = 1.
#' @param ssb.ss The spawning season (numeric). Default value = 1.
#' @param sr_model A character with the name of the model to simulate the recruitment process.
#' @param sr_params A named vector with the SR parameter values.
#' @param Fprop A vector with the same length as the number of seasons with the proportion of F by season. 
#'              By default the same proportion in all the seasons is assumed.
#' @param Fscan A vector with the F values to be simulated. 
#'              By default, will scan values between 0 and 4 (by 0.01 increments).
#' @param oldest The maximum age class to be considered. Default 100.
#' @param nrun The maximum number of years to project forward until equilibrium.
#' @param tol The desired accuracy.
#
#' @return A list with 2 elements runs and refpts
##'\itemize{
#'      \item{runs: is a data.frame containing a summary of the simulation results (year by year until equilibrium reached): 
#'          target F (ftgt), initial proportion of F in each season (fprop), expected biomass, B1+, ssb, recruitment, 
#'          spawning biomass per recruit (ssbpr), global fishing mortality (fbar), global F proportion by season (pF_s), 
#'          expected yields and their seasonal proportions by season (pY_s), yield per recruit (YPR), 
#'          harvest rates (hr_b1plus=catch/b1plus and hr_ssb=catch/ssb), covergence (TRUE or FALSE), 
#'          virgin status levels (rec0, ssb0, ssbpr0), status relative to virgin levels (ratioR, ratioB, ratioSPR), 
#'          slope (msy: slope = 0 at equilibrium) and rslope (f01: rslope = 0.1 at equilibrium).} 
#'      \item{refpts: is a data.frame with the different reference points estimates, giving equilibrium values for 
#'          some variables (ftgt, fprop_s1, fprop_s2, biomass, b1plus, ssb, rec, ssbpr, fbar).
#'          Showed reference levels are: msy, f0.1, %spr (spr.p), 
#'          and %B0 (F.pB0) for different percentages (p values) and F leading to 50%R0 (F.50R0).}
#'}
#'          
#'            
#' @seealso \code{\link{plotBRPsson}}   
#
#' @examples
#'
#' library(FLBEIA)
#' data(multistk) # object with 2 seasons and 3 iterations
#' 
#' stk <- trim( multistk, year=1, iter = 1) 
#' 
#' # loop for different catch proportions by season
#' for (p in seq(0.1,0.9,0.1)) {
#'   fruns <- brpsson( stk, B0=1e+05, R0=27489766, rec.ss=2, ssb.ss=2, 
#'                     sr_model="bevholt", sr_params=c( a = 29988835.109, b = 9090.909), 
#'                     Fprop = c(p, 1-p))
#' } 
#' 
#' plotBRPsson( fruns, pdfnm="stk_Fbar_vs_SPR.pdf")
#'  


#------------------------------------------------------------------------------#
# brpsson(obj) :: reference points with seasonal steps  
#------------------------------------------------------------------------------#

#' @rdname brpsson
brpsson <- function( stk, B0, R0, rec.ss=1, ssb.ss=1, sr_model, sr_params, 
                     Fprop = rep(1/dim(stk)[4], dim(stk)[4]), 
                     Fscan=seq(0,4,by=0.01), 
                     oldest=100, nrun = 200, tol = 1e-2) {
  
  require(dplyr)
  
  # Check
  
  if (ssb.ss>rec.ss) stop("'ssb.ss <= rec.ss' condition not met")
  
  
  # Dimensions
  
  nf <- length(Fscan)
  
  na <- dim(stk)[1]
  ns <- dim(stk)[4]
  ny <- nrun
  
  ages  <- as.numeric(dimnames(stk)[[1]])
  ages0 <- as.character(ages)
  ages1 <- as.character((ages[na]+1):(oldest+ages[1]-1))
  ages  <- c(ages0,ages1)
  if (ns > 1) { sson  <- as.numeric(dimnames(stk)[[4]]) } else sson <- 1
  year  <- 1:nrun
  
  fbar.range <- as.character(stk@range["minfbar"]:stk@range["maxfbar"])
  
  nage <- cage <- array(NA, c(oldest,ny,ns), dimnames = list(age=ages, year=year, season=sson))
  
  sage <- wagest <- wagec <- mage <- matage <- F.spwn <- M.spwn <-array(NA, c(oldest,ns), dimnames = list(age=ages, season=sson))
  
  ctot <- array(NA, c(ny,ns), dimnames = list(year=year, season=sson))
  
  SSB <- R <- structure(numeric(ny), names=year) * NA
  
  
  # Objects
  
  sage[ages0,]   <- stk@harvest[, drop = TRUE]                         # [na, ns]
  sage[ages1,]   <- rep(sage[na,],each=oldest-na)
  
  wagest[ages0,] <- apply(stk@stock.wt, c(1,4), mean)[drop = TRUE]     # [na, ns]
  wagest[ages1,] <- rep(wagest[na,],each=oldest-na)
  wagec[ages0,]  <- apply(stk@catch.wt, c(1,4), mean)[drop = TRUE]     # [na, ns]
  wagec[ages1,]  <- rep(wagec[na,],each=oldest-na)
  
  mage[ages0,]   <- apply(stk@m, c(1,4), mean)[drop = TRUE]            # [na, ns]
  mage[ages1,]   <- rep(mage[na,],each=oldest-na)
  matage[ages0,] <- apply(stk@mat, c(1,4), mean)[drop = TRUE]          # [na, ns]
  matage[ages1,] <- rep(matage[na,],each=oldest-na)
  
  
  F.spwn[ages0,] <- apply(stk@harvest.spwn, c(1,4), mean)[drop = TRUE] # [na, ns]
  F.spwn[ages1,] <- rep(F.spwn[na,],each=oldest-na)
  M.spwn[ages0,] <- apply(stk@m.spwn, c(1,4), mean)[drop = TRUE]       # [na, ns]
  M.spwn[ages1,] <- rep(M.spwn[na,],each=oldest-na)
  
  
  # Test different F values
  
  for (i in 1:nf) {
    
    ftgt <- Fscan[i] * Fprop
    
    fage <- t(ftgt * t(sage))
    
    zage.ssb <- fage[,ssb.ss] * F.spwn[,ssb.ss] + mage[,ssb.ss] * M.spwn[,ssb.ss] # [na]
    zage     <- fage + mage                                       # [na,ns]
    
    
    # 1st year
    
    nage[1,1,1:(rec.ss-1)] <- 0
    nage[1,1,rec.ss] <- R0                                             # - age 0
    if (rec.ss<ns) 
      for (ss in (rec.ss+1):ns) 
        nage[1,1,ss] <- nage[1,1,ss-1] * exp(-mage[1,ss-1])
    
    for (a in 2:length(ages))                                          # - intermediate ages
      for (ss in 1:ns)
        if (ss==1) {
          nage[a,1,ss] <- nage[a-1,1,ns] * exp(-mage[a-1,ns])
        } else {
          nage[a,1,ss] <- nage[a,1,ss-1] * exp(-mage[a,ss-1])
        }
    
    # # Alternative not to use a great number of age classes (i.e. older parameter)
    # nage[na,yr,1] <- nage[na-1,yr-1,ns] / (1- exp(-zage[na,ns]))
    
    # nage[length(ages),1,1] <- sum(nage[length(ages),1,1] *   # - plusgroup
    #                                      exp(-(0:(100-ages[length(ages)]))*sum(mage[length(ages),])))
    #
    # for (ss in 2:ns)
    #   nage[length(ages),1,ss] <- nage[length(ages),1,ss-1] * exp(-zage[length(ages),ss-1])
    
    
    SSB[1] <- sum( nage[,1,ssb.ss] * wagest[,ssb.ss] * matage[,ssb.ss] * exp(-zage.ssb))
    R[1]   <- R0
    
    cage[,1,] <- ctot[1,] <- 0
    
    
    # # Next seasons - 1st year
    # 
    # if (rec.ss<ns) 
    #   for (ss in (rec.ss+1):ns) 
    #     nage[,1,ss] <- nage[,1,ss-1] * exp(-zage[,ss-1])
    
    
    # Next years
    
    for (yr in 2:ny) {
      for (ss in 1:ns) {
        
        if (ss==1) {
          
          nage[-1,yr,ss] <- nage[-oldest,yr-1,ns] * exp(-zage[-oldest,ns])
          
          # nage[-c(1,na),yr,ss] <- nage[-c(na-1,na),yr-1,ns] * exp(-zage[-c(na-1,na),ns])
          # 
          # nage[na,yr,ss]       <- nage[na-1,yr-1,ns] * exp(-zage[na-1,ns]) + nage[na-1,yr-1,ns] * exp(-zage[na-1,ns])
          
        } else {
          
          nage[,yr,ss] <- nage[,yr,ss-1] * exp(-zage[,ss-1])
          
        }
        
        if (ss==rec.ss) {
          
          SSB[yr] <- sum( nage[,yr,ss] * wagest[,ss] * matage[,ss] * exp(-zage.ssb), na.rm = TRUE)
          
          # model call
          if(length(grep('~', sr_model)) == 0) { 
            model <- eval(call(sr_model))[[2]]
          } else # character but 'formula' 
            model <- formula(sr_model)
          
          datam <- as.list(sr_params)
          datam$ssb <- SSB[yr]
          
          R[yr]   <- c(eval(as.list(model)[[3]], datam)) 
          
          nage[1,yr,ss] <- R[yr]
          
        } else 
          nage[1,yr,ss] <- 0
      }
      
      cage[,yr,] <- fage/zage * nage[,yr,] * (1-exp(-zage))
      cage[is.nan(cage)] <- 0 
      
      ctot[yr,] <- apply(cage[,yr,]*wagec, 2, sum)
      
      if (abs(diff(R[c(yr-1,yr)]))<=tol) {
        conv <- TRUE
        break()
      }
      
    }
    
    
    # Summary output:
    
    biomass <- sum( nage[,yr,1] * wagest[,1]) 
    b1plus <- sum( nage[ages[ages!="0"],yr,1] * wagest[ages[ages!="0"],1])
    # ssb0   <- SSB[1]    # out %>% mutate(ssb0 = ssb[F==0])
    ssby   <- SSB[yr]
    # R0     <- R[1]      # out %>% mutate(rec0 = rec[F==0])
    Ry     <- R[yr]
    # ssbpr0 <- ssb0/r0   # out %>% mutate(ssbpr0 = ssbpr[F==0])
    ssbpr  <- ssby/Ry
    fbars  <- apply(fage[fbar.range,,drop=FALSE], 2, mean)
    fbary  <- sum(fbars)
    fbarp  <- fbars/sum(fbars)
    fbarp[is.nan(fbarp)] <- Fprop
    ctoty  <- sum(ctot[yr,])
    ctots  <- ctot[yr,]/ctoty
    ctots[is.nan(ctots)] <- 0
    ypr    <- ctoty/R[yr]
    hr_b1plus <- ctoty/b1plus
    hr_ssb    <- ctoty/ssby
    
    
    out.val <- c(Fscan[i], Fprop, biomass, b1plus, ssby, Ry, ssbpr, fbary, fbarp, ctoty, ctots, ypr, hr_b1plus, hr_ssb, conv)
    
    if (i==1) {
      out <- setNames(as.data.frame(t(out.val)), 
                      c("ftgt",paste("fprop",1:ns,sep="_s"),"biomass","b1plus","ssb","rec","ssbpr","fbar", paste("pF",1:ns,sep="_s"),"yield", paste("pY",1:ns,sep="_s"),
                        "YPR", "hr_b1plus", "hr_ssb", "conv"))
    } else
      out <- rbind( out, out.val)
    
    # out.val <- c( Fscan[i], b1plus, ssby, ssb0, ssby/ssb0, Ry, R0, Ry/R0, ssbpr, ssbpr0, ssbpr/ssbpr0, 
    #               fbary, fbarp, ctoty, ctots, ypr, hr_b1plus, hr_ssb, conv)
    # 
    # if (i==1) {
    #   
    #   out <- structure( out.val,
    #                     names=c("ftgt","b1plus","ssb","ssb0","ratSSB","rec","rec0","ratR","ssbpr","ssbpr0","pSPR",
    #                             "Fbar", paste("pF",1:ns,sep="_s"),"C", paste("pC",1:ns,sep="_s"), "YPR", "hr_b1plus", "hr_ssb", "conv"))
    # } else 
    #   out <- rbind(out, out.val)
    
    # stop if ratioB < 0.5% (i.e. = 0.005, that is when stock collapsed)
    if (ssby/SSB[1] < 0.005) break()
    
  } # end loop Fscan
  
  
  out <- out %>%
    mutate(rec0=rec[ftgt==0],
           ssb0=ssb[ftgt==0],
           ssbpr0=ssbpr[ftgt==0]) %>%
    mutate(ratioR=rec/rec0,
           ratioB=ssb/ssb0,
           ratioSPR=ssbpr/ssbpr0) # same as ssb/ssb0 if R is constant
  
  
  # Estimate slope and relative slope (Cy+1 - Cy)/(Fy+1/Fy)
  
  out <- out %>% arrange(ftgt) %>%
    mutate(cdiff = lead(yield, 1, order_by = ftgt) - yield, 
           fdiff = lead(fbar, 1, order_by = ftgt) - fbar, 
           slope = cdiff/fdiff, 
           rslope = slope/slope[ftgt==0]) %>%
    select(-cdiff, -fdiff)
  
  
  # Reference points
  
  aux <- out %>% group_by_at(vars(one_of(paste("fprop_s",1:ns,sep = "")))) #group_by(fprop_s1, fprop_s2)
  
  msy.val <- aux %>% filter(abs(slope - 0) == min(abs(slope - 0), na.rm = TRUE))
  f01.val <- aux %>% filter(abs(rslope - 0.1) == min(abs(rslope - 0.1), na.rm = TRUE))
  
  # If maximum yield not reached --> set MSY to NA
  if ( !( 0 >= min(aux$slope, na.rm = TRUE) & 0 <= max(aux$slope, na.rm = TRUE)) )
    msy.val <- msy.val %>%
                mutate_at(vars(-ftgt, -num_range("fprop_s",1:ns)), function(x) NA)
  
  refpts.val <- bind_rows( "virgin" = aux %>% filter(ftgt==0), 
                           "msy"    = msy.val, 
                           "f0.1"   = f01.val,
                           
                           "spr.20" = aux %>% filter(abs(ratioSPR - .20) == min(abs(ratioSPR - .20), na.rm = TRUE)),
                           "spr.30" = aux %>% filter(abs(ratioSPR - .30) == min(abs(ratioSPR - .30), na.rm = TRUE)), 
                           "spr.35" = aux %>% filter(abs(ratioSPR - .35) == min(abs(ratioSPR - .35), na.rm = TRUE)), 
                           "spr.40" = aux %>% filter(abs(ratioSPR - .40) == min(abs(ratioSPR - .40), na.rm = TRUE)), 
                           "spr.50" = aux %>% filter(abs(ratioSPR - .50) == min(abs(ratioSPR - .50), na.rm = TRUE)), 
                           
                           "F.20B0" = aux %>% filter(abs(ratioB - .20) == min(abs(ratioB - .20), na.rm = TRUE)), 
                           "F.30B0" = aux %>% filter(abs(ratioB - .30) == min(abs(ratioB - .30), na.rm = TRUE)), 
                           "F.35B0" = aux %>% filter(abs(ratioB - .35) == min(abs(ratioB - .35), na.rm = TRUE)), 
                           "F.40B0" = aux %>% filter(abs(ratioB - .40) == min(abs(ratioB - .40), na.rm = TRUE)), 
                           "F.50B0" = aux %>% filter(abs(ratioB - .50) == min(abs(ratioB - .50), na.rm = TRUE)), 
                           
                           "F.50R0" = aux %>% filter(abs(ratioR - .50) == min(abs(ratioR - .50), na.rm = TRUE)),
                           
                           .id = "refpt")
  
  msy09  <- refpts.val %>% filter(refpt=="msy") %>% .$yield * 0.90
  
  if (!is.na(msy09))
    f09msy <- aux %>% 
      filter(abs(yield - msy09) == min(abs(yield - msy09), na.rm = TRUE))
  else
    f09msy <- msy.val
  
  refpts.val <- bind_rows( refpts.val, "F.09msy" = f09msy  %>% mutate(refpt="F.90msy")) %>% 
    select(paste("fprop",1:ns,sep="_s"), refpt, fbar, yield, paste("pY",1:ns,sep="_s"), ssb, ratioB, ratioSPR, ratioR, hr_b1plus, hr_ssb)
  
  
  # keep which ftgt corresponds to fmsy and f01
  out <- out %>% mutate(fmsy = case_when(ftgt == msy.val$ftgt ~ TRUE,
                                         TRUE ~ FALSE), 
                         f01 = case_when(ftgt == f01.val$ftgt ~ TRUE,
                                         TRUE ~ FALSE))
  
  # refpts.val$refpt <- factor( refpts.val$refpt, 
  #                             levels = c("virgin", "msy", "f0.1", "spr.20", "spr.30", "spr.35", "spr.40", "spr.50", 
  #                                        "F.20B0", "F.30B0", "F.35B0", "F.40B0", "F.50B0", "F.50R0", "F.90msy"))
  
  
  # Return output
  
  return( list( runs = as.data.frame(out), refpts = refpts.val))
  
}


#------------------------------------------------------------------------------#
# plotBRPsson (obj) :: YPR plot from brpsson output  
#------------------------------------------------------------------------------#


#' @rdname plotBRPsson
#' 
#' @title YPR plots for a seasonal model
#' 
#' @description Function for generating YPR plots for reference points when having a seasonal model
#'
#' @param obj The output of brpsson function or a data.frame with the same structure as the runs object in the 
#'            output list of brpsson function.
#' @param pdfnm The name of the pdf document where plots are going to be saved.
#
#' @return A pdf for each stock with plots.
#'          
#' @seealso \code{\link{brpsson}}          
#
# @examples
# 
# 

plotBRPsson <- function( obj, pdfnm="Fbar_vs_SPR.pdf") {
  
  
  require(ggplot2)
  
  
  # Check object (can be brpsson output or data.frame)
  
  if (inherits(obj, "list")) obj <- obj$runs
  
  
  # Check object has correct colunms
  
  ind <- c("ftgt", "fprop_s1", "fbar", "yield", "pY_s1", "YPR", "ratioSPR") 
  
  imiss <- !ind %in% names(obj)
  
  if (any(imiss))
    stop(paste("Variables ", paste("'",ind[imiss],"'", collapse=", ", sep=""), " are missing.", sep=""))
  
  
  pdf(pdfnm, onefile=T)
  
  # Fbar vs. %SPR & catch
  secaxis.rat <- 1e+5/1.5
  
  p <- ggplot(obj, aes(x = fbar, y = ratioSPR)) + 
    geom_line(aes(y = ratioSPR, color="ratioSPR")) + 
    geom_line(aes(y = yield*1.5e-5, color = "yield")) + 
    facet_wrap( . ~ fprop_s1, ncol = 3) + 
    scale_y_continuous(sec.axis = sec_axis(~.*secaxis.rat, name = "yield")) + 
    scale_color_manual(name = "", values = c("ratioSPR"="black", "yield"="blue")) + 
    labs(y = "%SPR",
         x = "Fbar") + 
    theme(text=element_text(size=12),
          strip.text=element_text(size=12),
          title=element_text(size=16,face="bold"),
          axis.title.y.right=element_text(color="blue"), 
          legend.direction = "horizontal",
          legend.position = "top", 
          legend.title = element_blank()) +
    ggtitle("YPR (by fprop_s1)")
  
  # Add F0.1
  p <- p + geom_segment( aes(x = ftgt, y = 0, xend = ftgt, yend = yield/secaxis.rat), linetype=3, colour="red", 
                         data = obj %>% filter(f01 == TRUE)) +
    geom_segment( aes(x = 0, y = yield/secaxis.rat, xend = ftgt, yend = yield/secaxis.rat), linetype=3, colour="red", 
                  data = obj %>% filter(f01 == TRUE)) + 
    geom_text( aes( x=ftgt-0.3, y=0.05), label="F0.1", size = 2, colour="red", 
               data = obj %>% filter(f01 == TRUE))
  
  # Add Fmsy
  p <- p + geom_segment( aes(x = ftgt, y = 0, xend = ftgt, yend = yield/secaxis.rat), linetype=3, colour="red", 
                         data = obj %>% filter(fmsy == TRUE)) +
    geom_segment( aes(x = 0, y = yield/secaxis.rat, xend = ftgt, yend = yield/secaxis.rat), linetype=3, colour="red", 
                  data = obj %>% filter(fmsy == TRUE)) + 
    geom_text( aes( x=ftgt+0.3, y=0.05), label="Fmsy", size = 2, colour="red", 
               data = obj %>% filter(fmsy == TRUE))
  
  print(p)
  
  
  # Fbar vs. %SPR & catch & pC_s1
  
  p <- ggplot(obj, aes(x = fbar, y = ratioSPR)) + 
    geom_line(aes(y = ratioSPR, color="ratioSPR")) + 
    geom_line(aes(y = yield*1.5e-5, color = "yield")) + 
    geom_line(aes(y = pY_s1, color = "pY_s1")) +
    facet_wrap( . ~ fprop_s1, ncol = 3) + 
    scale_y_continuous(sec.axis = sec_axis(~.*secaxis.rat, name = "yield")) + 
    scale_color_manual(name = "", values = c("ratioSPR"="black", "yield"="blue", "pY_s1"="green"), 
                       labels = c("ratioSPR"="%SPR", "yield"="yield", "pY_s1"="%yield_sem1")) + 
    labs(y = "%",
         x = "Fbar (1-3)") + 
    theme(text=element_text(size=12),
          strip.text=element_text(size=12),
          title=element_text(size=16,face="bold"),
          axis.title.y.right=element_text(color="blue"), 
          legend.direction = "horizontal", 
          legend.position = "top") +
    theme_bw() +
    ggtitle("YPR (by pF_sem1)")
  
  print(p)
  
  dev.off()
  
  
}
