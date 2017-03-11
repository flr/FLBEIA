#-------------------------------------------------------------------------------
#
#' Plots with biological data 
#' 
#' For each stock, return a pdf with plots using FLBiols object.
#'
#' @details
#'\itemize{
#'      \item{Each pdf contains biomass in numbers at age, mean weight at age, fecundity, natural mortality, maturity, spawning, recruitment and spawning stock biomass} 
#'}


#' @param biols A FLBiols object 
#' @param pdfnm The name for the pdf document will be stock's name and pdfnm separated by a line.
#
#' @return A pdf for each stock with plots.

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(ggplot2)
#' data(one)
#' s0 <- FLBEIA(biols = oneBio,       # FLBiols object with one FLBiol element for stk1.
#'                SRs = oneSR,        # A list with one FLSRSim object for stk1.
#'                BDs = NULL,         # No Biomass Dynamic populations in this case.
#'             fleets = oneFl,        # FLFleets object with on fleet.
#'             covars = NULL,         # covars not used
#'            indices = NULL,         # indices not used 
#'             advice = oneAdv,       # A list with two elements 'TAC' and 'quota.share'
#'          main.ctrl = oneMainC,     # A list with one element to define the start and end of the simulation.
#'         biols.ctrl = oneBioC,      # A list with one element to select the model to simulate the stock dynamics.
#'        fleets.ctrl = oneFlC,       # A list with several elements to select fleet dynamic models and store additional parameters.
#'        covars.ctrl = NULL,         # covars control not used 
#'           obs.ctrl = oneObsC,      # A list with one element to define how the stock observed ("PerfectObs").
#'        assess.ctrl = oneAssC,      # A list with one element to define how the stock assessment model used ("NoAssessment").
#'        advice.ctrl = oneAdvC) 
#' plotFLBiols(s0$biols, 's0')
#' }


###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables in a biols object
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................
plotFLBiols <- function(biols,pdfnm){
    
  names.biols <- names(biols)
  path.pdf <- ''
  
  for(i in 1:length(names.biols)){
    
    pdf(paste(path.pdf,names.biols[i],'-',pdfnm,'.pdf',sep=''))
    
    biol <- biols[[i]]
 

    biol.n.df <- as.data.frame(biol@n)
    biol.n.df$variable <- 'n'
    biol.n.df$indicator <- names.biols[i] 
    biol.n.df$age <- factor(biol.n.df$age)
    
    biol.wt.df <- as.data.frame(biol@wt)
    biol.wt.df$variable <- 'wt'
    biol.wt.df$indicator <- names.biols[i] 
    biol.wt.df$age <- factor(biol.wt.df$age)
    
    df <- rbind(biol.n.df,biol.wt.df) 
    df$stock <- names.biols[i] 
    temp <- aggregate(data ~ age+year+stock+indicator+variable, data = df, 
                      mean, na.rm=TRUE,na.action="na.pass")
    
    p <- ggplot(data=temp, aes(x=year, y=data, fill=age))  + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(variable~indicator,scales=c("free_y"))
    print(p)   
    
    df <- NULL
    biol.m.df <- as.data.frame(biol@m)
    biol.m.df$variable <- 'm'
    biol.m.df$indicator <- names.biols[i] 
    biol.m.df$age <- factor(biol.n.df$age)
    
    biol.fec.df <- as.data.frame(fec(biol))
    biol.fec.df$variable <- 'fec'
    biol.fec.df$indicator <- names.biols[i] 
    biol.fec.df$age <- factor(biol.fec.df$age)

    biol.mat.df <- as.data.frame(mat(biol))
    biol.mat.df$variable <- 'mat'
    biol.mat.df$indicator <- names.biols[i] 
    biol.mat.df$age <- factor(biol.mat.df$age)
    
            
    biol.spwn.df <- as.data.frame(spwn(biol))
    biol.spwn.df$variable <- 'spwn'
    biol.spwn.df$indicator <- names.biols[i] 
    biol.spwn.df$age <- factor(biol.spwn.df$age)
    
    df <- rbind(biol.m.df,biol.fec.df,biol.mat.df,biol.spwn.df) 
    df$stock <- names.biols[i] 
    temp <- aggregate(data ~ age + year+variable+stock+indicator, data = df, 
                      mean, na.rm=TRUE,na.action="na.pass")
    
    p <- ggplot(data=temp, aes(x=year, y=data, fill=age))  + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(variable~indicator,scales=c("free_y"))
    print(p)  
    
    df<- NULL
    biol.ssb.df <- as.data.frame(ssb(biol))
    biol.ssb.df$variable <- 'ssb'
    biol.ssb.df$indicator <- names.biols[i]
    biol.ssb.df$age <- factor(biol.ssb.df$age)
    rec <- biol@n[1,,1,1,,]
    n.ss <- dim(biol@n)[4]
    n.unit <- dim(biol@n)[3]
    if(n.ss>2 & n.unit>2){
      for(k in 2:length(n.ss)){
        rec <- rec+biol@n[1,,k,k,,]}}
    biol.rec.df <- as.data.frame(rec)
    biol.rec.df$variable <- 'rec'
    biol.rec.df$indicator <- names.biols[i]
    biol.rec.df$age <- factor(biol.rec.df$age)
    
    df <- rbind(biol.ssb.df,biol.rec.df) 
    df$stock <- names.biols[i] 
    temp <- aggregate(data ~ age + year+variable+stock+indicator, data = df, 
                      mean, na.rm=TRUE,na.action="na.pass")
    
    p <- ggplot(data=temp, aes(x=year, y=data, fill=age))  + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(variable~indicator,scales=c("free_y"))
    print(p)  
    
    
    dev.off()   
    
  }
}
    


