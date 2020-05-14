#' Plots with fleets data 
#' 
#' Return a pdf with plots using FLBEIA object (FLFleets and covars).
#'
#' @param obj An FLBEIA object. 
#' @param probs A numeric vector with the probabilities used to calculate the quantiles.
#' @param pdfnm The name for the pdf document will be "Eco" and pdfnm separated by a line.
#
#' @return A pdf with capacity, costs, effort, profits by fleet.

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(ggplot2)
#' 
#' data(res_flbeia)
#' plotEco(oneRes, pdfnm='one')
#' }


###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables from the data frame coming from ecoSum function
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................

plotEco <- function(obj,probs = c(0.95,0.5,0.05),pdfnm="bc"){

 names.fl <- names(obj$fleets)

 eco <- fltSum(obj, flnms = names.fl, years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = FALSE)
   
 

  path.pdf <- ''
   pdf(paste(pdfnm,'_',path.pdf,'Eco','.pdf',sep='')) 
  for(i in 1:length(names.fl)){
    

    fleet <- obj$fleets[[i]]
    
      #capacity   

    
    res <- aggregate(capacity ~ year+ fleet, eco, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    res$fleet <- as.factor(res$fleet)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50,fill=.data$fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=fleet), alpha=0.3) + 
      ggtitle("Capacity")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
    #costs
    res <- aggregate(costs ~ year+ fleet, eco, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=fleet), alpha=0.3) + 
      ggtitle("Costs")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)

    #effort
    res <- aggregate(effort ~ year+ fleet, eco, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=fleet), alpha=0.3) + 
      ggtitle("Effort")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
    #grossSurplus
    res <- aggregate(grossSurplus ~ year+ fleet, eco, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=fleet), alpha=0.3) + 
      ggtitle("Gross-Surplus")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
  }  
 dev.off()
}


