###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables in a biols object
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................
plot.biols <- function(biols,pdfnm){
  
  require(ggplot2)
  require(plyr)
  require(FLCore)
  
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
    
    biol.fec.df <- as.data.frame(biol@fec)
    biol.fec.df$variable <- 'fec'
    biol.fec.df$indicator <- names.biols[i] 
    biol.fec.df$age <- factor(biol.fec.df$age)
        
    biol.spwn.df <- as.data.frame(biol@spwn)
    biol.spwn.df$variable <- 'spwn'
    biol.spwn.df$indicator <- names.biols[i] 
    biol.spwn.df$age <- factor(biol.spwn.df$age)
    
    df <- rbind(biol.m.df,biol.fec.df,biol.spwn.df) 
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
    
    biol.rec.df <- as.data.frame(rec(biol))
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
    


