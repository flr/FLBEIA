#' Plots with fltStkSum data 
#' 
#' Return a pdf with plots with the outputs fltStkSum, usint the output of FLBEIA.
#'
#' @param obj FLBEIA output
#' @param pdfnm The name for the pdf document will be 'fltStkSum-' and pdfnm.
#
#' @return One pdf with plots on landings, discards,catch,price,quotaUpt,tacshare,discRat,quota by fleet and stock.

#' @examples
#'\dontrun{
#' library(FLBEIA)
#' library(ggplot2)
#' data(res_flbeia)
#' plotfltStkSum(obj=oneRes,pdfnm = "oneRes")
#' }
#' 
#.......................................................
#....................FUNCTIONS..........................

plotfltStkSum <- function(obj,pdfnm){
  
 fltStkSum.data  <-fltStkSum(obj, flnms = "all", years = dimnames(obj$biols[[1]]@n)$year,
                           byyear = TRUE, long = TRUE, scenario = "bc")

  path.pdf <- ''
  pdf(paste(pdfnm,'_',path.pdf,'fltStkSum','.pdf',sep=''))
  
  indicator <- unique(fltStkSumQ(fltStkSum.data)$indicator)
  
  for(i in 1:length(indicator)){
    
    data <- fltStkSumQ(fltStkSum.data)
    sub.data <- data[data$indicator==indicator[i],]

    # temp <- aggregate(value ~ year+ fleet+stock, 
    #                   data = sub.data, mean , na.rm=TRUE,na.action="na.pass") 
    # temp$year <- as.numeric(as.character(temp$year))
    sub.data$year <- as.numeric(sub.data$year)
    p <- ggplot( data=sub.data, aes(x=year, y=q50, color=stock)) + 
      geom_line() + theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=year, ymin=q05, ymax=q95, fill=stock), alpha=0.5) + 
      facet_wrap(~ fleet, scales="free") + 
      ylab(indicator[i])
    print(p)
    

  }  
 dev.off()
}


