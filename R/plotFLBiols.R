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
#' @param prob a numeric vector with the probabilities used to calculate the quantiles. 
#' @return A pdf for each stock with plots.

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(ggplot2)
#' data(res_flbeia)
#' plotFLBiols(oneRes$biols, pdfname='oneRes')
#' }


###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables in a biols object
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................
plotFLBiols <- function(biols,prob = c(0.95,0.5,0.05),pdfnm="bc"){
    
  names.biols <- names(biols)
  path.pdf <- ''
  
  for(i in 1:length(names.biols)){
    
    pdf(paste(pdfnm,'_',path.pdf,names.biols[i],'.pdf',sep=''))
    
    biol <- biols[[i]]
    
    biol.sl.df <- as.data.frame(biol@n)
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
       ggtitle("n")+theme(plot.title = element_text(hjust = 0.5)) 
      
    print(p)
    
   
    biol.sl.df <- as.data.frame(biol@wt)
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("wt")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(biol@m)
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("m")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(fec(biol))
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("fec")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(mat(biol))
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("mat")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(spwn(biol))
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("spwn")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
 
    biol.sl.df <- as.data.frame(ssb(biol))
    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("ssb")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    rec <- biol@n[1,,1,1,,]
    n.ss <- dim(biol@n)[4]
    n.unit <- dim(biol@n)[3]
    if(n.ss>2 & n.unit>2){
      for(k in 2:length(n.ss)){
        rec <- rec+biol@n[1,,k,k,,]}}
    biol.sl.df <- as.data.frame(rec)

    res <- aggregate(data ~ year + season+age, biol.sl.df, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[4:(4+length(prob)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=year, y=q50, fill=age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=age), alpha=0.3) + 
      ggtitle("rec")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)

    dev.off()   
    
  }
}
    


