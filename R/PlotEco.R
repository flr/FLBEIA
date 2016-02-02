###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables from the data frame coming from ecoSum function
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................

plot.eco <- function(fleets,pdfnm){
  
 require(ggplot2)
 require(plyr)
 require(FLCore)

 names.fl <- names(fleets)
 eco  <- ecoSum(fleets, flnms= names.fl, year=dimnames(fleets[[1]]@effort)$year)
 

  path.pdf <- ''
   pdf(paste(path.pdf,'Eco','-',pdfnm,'.pdf',sep='')) 
  for(i in 1:length(names.fl)){
    

    fleet <- fleets[[i]]
    
      #capacity   
    temp <- aggregate(capacity ~ year+ fleet, 
                      data = eco, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=capacity)) + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(.~fleet,scales=c("free_y"))     

    print(p)
    
    #costs
    temp <- aggregate(costs ~ year+fleet, 
                      data = eco, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=costs)) + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(.~fleet,scales=c("free_y"))     
    
    print(p)

    #effort
    temp <- aggregate(effort ~ year+ fleet, 
                      data = eco, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=effort)) + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(.~fleet,scales=c("free_y"))     
    
    print(p)
    
    #profits
    temp <- aggregate(profits~ year+ fleet, 
                      data = eco, mean , na.rm=TRUE,na.action="na.pass") 
    temp$year <- as.numeric(as.character(temp$year))
    
    p <- ggplot(data=temp, aes(x=year, y=profits)) + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(.~fleet,scales=c("free_y"))     
    
    print(p)

  }  
 dev.off()
}


