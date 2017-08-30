#' Plots with fleets data 
#' 
#' Return a pdf with plots using FLBEIA object (FLFleets and covars).
#'
#' @param obj An FLBEIA object. 
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

plotEco <- function(obj,prob = c(0.95,0.5,0.05),pdfnm="bc"){

 names.fl <- names(obj$fleets)
 # eco  <- ecoSum_damara(obj$fleets, flnms= names.fl, years=dimnames(obj$fleets[[1]]@effort)$year)
 eco <- fltSum(obj, flnms = names.fl, years = dimnames(obj$biols[[1]]@n)$year, byyear = TRUE, long = FALSE)
   
 

  path.pdf <- ''
   pdf(paste(pdfnm,'_',path.pdf,'Eco','.pdf',sep='')) 
  for(i in 1:length(names.fl)){
    

    fleet <- obj$fleets[[i]]
    
      #capacity   

    
    res <- aggregate(capacity ~ year+ fleet, eco, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[3:(3+length(prob)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    res$fleet <- as.factor(res$fleet)
    p <- ggplot( data=res, aes(x=year, y=q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=fleet), alpha=0.3) + 
      ggtitle("capacity")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
    #costs
    res <- aggregate(costs ~ year+ fleet, eco, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[3:(3+length(prob)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    p <- ggplot( data=res, aes(x=year, y=q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=fleet), alpha=0.3) + 
      ggtitle("costs")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)

    #effort
    res <- aggregate(effort ~ year+ fleet, eco, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[3:(3+length(prob)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    p <- ggplot( data=res, aes(x=year, y=q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=fleet), alpha=0.3) + 
      ggtitle("effort")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
    #profits
    res <- aggregate(profits ~ year+ fleet, eco, quantile, prob = prob, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(prob,3, nchar(prob)))==1, paste(substr(prob,3, nchar(prob)), 0, sep = ""), substr(prob,3, nchar(prob))), sep = "")
    names(res)[3:(3+length(prob)-1)] <- nms
    res$year <- as.numeric(as.character(res$year))
    p <- ggplot( data=res, aes(x=year, y=q50,fill=fleet)) + 
      geom_line() +geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=fleet), alpha=0.3) + 
      ggtitle("profits")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    
  }  
 dev.off()
}


