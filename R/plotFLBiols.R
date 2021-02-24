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
#' @param probs a numeric vector with the probabilities used to calculate the quantiles. 
#' @return A pdf for each stock with plots.

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(ggplot2)
#' data(res_flbeia)
#' plotFLBiols(oneRes$biols, pdfnm='oneRes')
#' }


###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables in a biols object
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................
plotFLBiols <- function(biols,probs = c(0.95,0.5,0.05),pdfnm="bc",u=1,ss=1){
    
  names.biols <- names(biols)
  path.pdf <- ''
  
  for(i in 1:length(names.biols)){
    
    pdf(paste(pdfnm,'_',path.pdf,names.biols[i],'.pdf',sep=''))
    
    biol <- biols[[i]]
    
    biol.sl.df <- as.data.frame(biol@n[,,u,ss])
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
       ggtitle(paste("n (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
      
    print(p)
    
    biol.sl.df <- as.data.frame(unitSums(biol@n[,,,ss]))
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("n (season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
   
    biol.sl.df <- as.data.frame(biol@wt[,,u,ss])
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("wt (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(biol@m[,,u,ss])
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("m (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(fec(biol)[,,u,ss])
    res <- aggregate(data ~ year +age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("fec (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(mat(biol)[,,u,ss])
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("mat (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(spwn(biol)[,,u,ss])
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("spwn (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
 
    biol.sl.df <- as.data.frame(ssb(biol)[,,u,ss])
    res <- aggregate(data ~ year +age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("ssb (unit",u,"& season",ss,")"))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    biol.sl.df <- as.data.frame(unitSums(ssb(biol)[,,,1]))
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle("Total ssb (beggining of the year) ")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    rec <- biol@n[1,,u,ss,,]

    biol.sl.df <- as.data.frame(rec)

    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle("rec")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
    rec <-seasonSums(unitSums(biol@n[1,,,,,]))
    
    biol.sl.df <- as.data.frame(rec)
    
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle("Total rec per year")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)

    dev.off()   
    
  }
}
    


