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
#' 
#' plotFLBiols(multiRes$biols, pdfnm='multiRes', season = 2) # only specific season
#' }


###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia, Sonia Sanchez-Maroño 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables in a biols object
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................
plotFLBiols <- function(biols, probs = c(0.95,0.5,0.05), pdfnm = "bc", season = "all"){
    
  if(length(season)>1) stop("'season' argument must be 'all' or specific season")
  
  names.biols <- names(biols)
  path.pdf <- ''
  
  # specific season naming
  ss.nam    <- ifelse(season != "all", paste0("_ss",season), "")
  ss.legend <- paste0("(season: ",season,")")
  
  if (season == "all") { ss <- 1:dims(biols[[1]])$season } else ss <- season
  
  for(i in 1:length(names.biols)){
    
    pdf(paste(pdfnm,'_',path.pdf,names.biols[i],ss.nam,'.pdf',sep=''))
    
    biol <- biols[[i]]
    
    # Numbers-at-age by unit
    
    biol.sl.df <- as.data.frame(biol@n[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
       ggtitle(paste("n", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    # Numbers-at-age all units
    
    biol.sl.df <- as.data.frame(unitSums(biol@n[,,,ss]))
    res <- aggregate(data ~ year + age + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:3], data.frame(res[,4]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[4:(4+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      facet_wrap(. ~ season) +
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("n", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    # Weights-at-age by unit
    
    biol.sl.df <- as.data.frame(biol@wt[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("wt", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    # Natural mortality-at-age by unit
    
    biol.sl.df <- as.data.frame(biol@m[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("m", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    # Fecundity-at-age by unit
    
    biol.sl.df <- as.data.frame(fec(biol)[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("fec", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    # Maturity-at-age by unit
    
    biol.sl.df <- as.data.frame(mat(biol)[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("mat", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    # Percentage of mortality before spawning-at-age by unit
    
    biol.sl.df <- as.data.frame(spwn(biol)[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("spwn", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
 
    # SSB by unit
    
    biol.sl.df <- as.data.frame(ssb(biol)[,,,ss])
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("ssb", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    # SSB at the beginning of the year
    
    biol.sl.df <- as.data.frame(unitSums(ssb(biol)[,,,1]))
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle("Total ssb (beggining of the year) ")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)
    
    # Recruitment by unit
    
    rec <- biol@n[1,,,ss,,]

    biol.sl.df <- as.data.frame(rec)
    res <- aggregate(data ~ year + age + unit + season, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:4], data.frame(res[,5]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[5:(5+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    res$unit   <- paste0("u = ", res$unit)
    res$season <- paste0("ss = ", res$season)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle(paste("rec", ss.legend))+theme(plot.title = element_text(hjust = 0.5)) 
    
    if (length(ss) > 1) {p <- p + facet_grid(season ~ unit)} else p <- p + facet_wrap(. ~ unit)
    
    print(p)
    
    
    # Total recrutiment
    
    rec <- seasonSums(unitSums(biol@n[1,,,,,]))
    
    biol.sl.df <- as.data.frame(rec)
    
    res <- aggregate(data ~ year + age, biol.sl.df, quantile, probs = probs, na.rm=T)
    res <- cbind(res[,1:2], data.frame(res[,3]))
    nms <- paste('q',ifelse(nchar(substr(probs,3, nchar(probs)))==1, paste(substr(probs,3, nchar(probs)), 0, sep = ""), substr(probs,3, nchar(probs))), sep = "")
    names(res)[3:(3+length(probs)-1)] <- nms
    res$age <- as.factor(res$age)
    p <- ggplot( data=res, aes(x=.data$year, y=.data$q50, fill=.data$age)) + 
      geom_line() + geom_point(size=2,shape=21)+ theme_bw() + 
      theme(text=element_text(size=10),
            title=element_text(size=10,face="bold"),
            strip.text=element_text(size=10)) + 
      geom_ribbon(aes(x=.data$year, ymin=.data$q05, ymax=.data$q95, fill=.data$age), alpha=0.3) + 
      ggtitle("Total rec per year")+theme(plot.title = element_text(hjust = 0.5)) 
    
    print(p)

    dev.off()   
    
  }
}
    


