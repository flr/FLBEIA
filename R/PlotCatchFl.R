###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables from the data frame coming from catchFlSum function
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................

plot.catchFl <- function(fleets,advice,pdfnm){
  
 require(ggplot2)
 require(plyr)
 require(FLCore)

 names.fl <- names(fleets)
 catchFl  <- catchFlSum(fleets,advice,flnms= names.fl,
                           stknms= 'all', years=dimnames(fleets[[1]]@effort)$year)
 

  path.pdf <- ''
  pdf(paste(path.pdf,'Catch','-',pdfnm,'.pdf',sep=''))
  for(i in 1:length(names.fl)){
    

    fleet <- fleets[[i]]
    
      #landings    
    temp <- aggregate(landings ~ year+ fleet+stock, 
                      data = catchFl, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=landings, fill=stock)) + 
      geom_area(colour="black", size=.2, alpha=.4)+ 
      facet_grid(.~fleet,scales=c("free_y"))     

    print(p)
    
    #discards
    temp <- aggregate(discards ~ year+ fleet+stock, 
                      data = catchFl, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=discards, fill=stock)) + 
      geom_area(colour="black", size=.2, alpha=.4)+ 
      facet_grid(.~fleet,scales=c("free_y"))     
    
    print(p)

    #price
    temp <- aggregate(price ~ year+ fleet+stock, 
                      data = catchFl, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=price, fill=stock))  + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(.~fleet,scales=c("free_y"))     
    
    print(p)
    

  }  
 dev.off()
}


