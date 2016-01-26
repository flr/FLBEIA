###############################################################################
# AUTHOR(DATE):       Agurtzane Urtizberea, Dorleta Garcia 
# RESEARCH INSTITUTE: AZTI-TECNALIA                      
# TITLE:        Check variables in a fleets object
# NOTE #1:      Return plots
###############################################################################
#.......................................................
#....................FUNCTIONS..........................
#total.catch.stock.df function returns a data frame with
#the total catch per stock for a fleet object.
#....................................................
#total.landings.stock.df function returns a data frame with
#the total landings per stock for a fleet object
#....................................................
#total.discards.stock.df function returns a data frame with
#the total discards per stock for a fleet object
#....................................................
#plot.fleets makes the plots of each variable in fleet



total.catch.stock.df <- function(fleet){

  nms.metiers <- names(fleet@metiers)
  nms.stock <- NULL
  for(i in 1:length(nms.metiers)){
    nms.stock <- c(nms.stock,names(fleet@metiers[[i]]@catches))
  }
  nms.stock <- unique(nms.stock)
  df.catch <- NULL
  
  for(j in 1:length(nms.stock)){
    nm.stock <- nms.stock[j]
    #catch.stock <- fleet@metiers[[1]]@catches[[nm.stock]@landings
    catch.stock <- 0
    for(i in 1:length(nms.metiers)){
      nm.met.stock <- names(fleet@metiers[[i]]@catches)
      if(nm.stock %in% nm.met.stock){
        catch.stock <- catch.stock +
          fleet@metiers[[i]]@catches[[nm.stock]]@landings +
          fleet@metiers[[i]]@catches[[nm.stock]]@discards
      }else{
        catch.stock <- catch.stock     
      }
    }
   #   list.catch[[nms.stock[j]]] <- catch.stock 
    df <- as.data.frame(catch.stock)
    df$species <- rep(nms.stock[j],dim(df)[1])
    df.catch <- rbind(df.catch,df)    
  }
return(df.catch)
}

total.landings.stock.df <- function(fleet){
  
  nms.metiers <- names(fleet@metiers)
  nms.stock <- NULL
  for(i in 1:length(nms.metiers)){
    nms.stock <- c(nms.stock,names(fleet@metiers[[i]]@catches))
  }
  nms.stock <- unique(nms.stock)
  df.landings <- NULL
  
  for(j in 1:length(nms.stock)){
    nm.stock <- nms.stock[j]
    #catch.stock <- fleet@metiers[[1]]@catches[[nm.stock]@landings
    landings.stock <- 0
    for(i in 1:length(nms.metiers)){
      nm.met.stock <- names(fleet@metiers[[i]]@catches)
      if(nm.stock %in% nm.met.stock){
        landings.stock <- landings.stock +
          fleet@metiers[[i]]@catches[[nm.stock]]@landings 
      }else{
        landings.stock <- landings.stock     
      }
    }
    #   list.catch[[nms.stock[j]]] <- catch.stock 
    df <- as.data.frame(landings.stock)
    df$species <- rep(nms.stock[j],dim(df)[1])
    df.landings<- rbind(df.landings,df)    
  }
  return(df.landings)
}



total.discards.stock.df <- function(fleet){
  
  nms.metiers <- names(fleet@metiers)
  nms.stock <- NULL
  for(i in 1:length(nms.metiers)){
    nms.stock <- c(nms.stock,names(fleet@metiers[[i]]@catches))
  }
  nms.stock <- unique(nms.stock)
  df.discards <- NULL
  
  for(j in 1:length(nms.stock)){
    nm.stock <- nms.stock[j]
    #catch.stock <- fleet@metiers[[1]]@catches[[nm.stock]@landings
    discards.stock <- 0
    for(i in 1:length(nms.metiers)){
      nm.met.stock <- names(fleet@metiers[[i]]@catches)
      if(nm.stock %in% nm.met.stock){
        discards.stock <- discards.stock +
          fleet@metiers[[i]]@catches[[nm.stock]]@discards
      }else{
        discards.stock <- discards.stock     
      }
    }
    #   list.catch[[nms.stock[j]]] <- catch.stock 
    df <- as.data.frame(discards.stock)
    df$species <- rep(nms.stock[j],dim(df)[1])
    df.discards<- rbind(df.discards,df)    
  }
  return(df.discards)
}
#........................................................
#........................................................



plot.fleets <- function(fleets,pdfnm){
  
 require(ggplot2)
 require(plyr)
 require(FLCore)

  names.fl <- names(fleets)
  path.pdf <- ''
  
  for(i in 1:length(names.fl)){
    
    pdf(paste(path.pdf,names.fl[i],'-',pdfnm,'.pdf',sep=''))
  
    fleet <- fleets[[i]]
    
    #TOTAL CATCH, LANDINGS AND DISCARDS
    total.catch.df <- total.catch.stock.df(fleet)
    total.catch.df$indicator <- 'catch'
        
    total.landings.df <- total.landings.stock.df(fleet)
    total.landings.df$indicator <- 'landings'
    
    total.discards.df <- total.discards.stock.df(fleet)
    total.discards.df$indicator <-'discards'
    
    df <- rbind(total.catch.df,total.landings.df,total.discards.df)
    df$fleet <- rep(paste(names.fl[i],dim(df)[1]))
    
    temp <- aggregate(data ~ age + year+species+indicator+fleet, data = df, mean)
    temp$indicator <- factor(temp$indicator)
    temp$species <- factor(temp$species)
    temp$fleet <- factor(temp$fleet)   
    p <- ggplot(data=temp, aes(x=year, y=data, fill=species))  + 
      geom_area(colour="black", size=.2, alpha=.4)+ 
      facet_grid(indicator~fleet,scales=c("free_y"))
    print(p)
    
    #EFFORT, FCOST,CAPACITY,CREWSHARE
    df <- NULL
    effort.df <- as.data.frame(fleet@effort)
    effort.df$variable <- 'effort'
    effort.df$indicator <- names.fl[i] 
   
    fcost.df <- as.data.frame(fleet@fcost)
    fcost.df$variable <- 'fcost'
    fcost.df$indicator <- names.fl[i] 
    
    capacity.df <- as.data.frame(fleet@capacity)
    capacity.df$variable <- 'capacity'
    capacity.df$indicator <- names.fl[i] 
    
    crewshare.df <- as.data.frame(fleet@crewshare)
    crewshare.df$variable <- 'crewshare'
    crewshare.df$indicator <- names.fl[i] 

    df <- rbind(df,effort.df,fcost.df,capacity.df,crewshare.df) 
 
    temp <- aggregate(data ~ age + year+area+variable+indicator, data = df, mean)
    
    p <- ggplot(data=temp, aes(x=year, y=data, fill=indicator))  + geom_line() +
      geom_point(size=2, shape=21)+
      facet_grid(variable~indicator,scales=c("free_y"))
    print(p)          
    
    #PER METIER
    nms.metiers <- names(fleet@metiers)
    
    df <- NULL
    for(k in 1:length(nms.metiers)){
      
      #EFFSHARE, VCOST
      effshare.df <- as.data.frame(fleet@metiers[[k]]@effshare)
      effshare.df$metier <- nms.metiers[k]
      effshare.df$indicator <- 'effshare'
            
      vcost.df <- as.data.frame(fleet@metiers[[k]]@effshare)
      vcost.df$metier <- nms.metiers[k]
      vcost.df$indicator <- 'vcost'
      df <- rbind(df,effshare.df,vcost.df) 
      
    }
      df$fleet <- names.fl[i]
    temp <- aggregate(data ~ age + year+area+fleet+metier+indicator, data = df, mean)
    
      p <- ggplot(data=temp, aes(x=year, y=data, fill=metier))  + geom_line() +
        geom_point(size=2, shape=21)+
        facet_grid(indicator~fleet,scales=c("free_y"))
      print(p)          

    nms.stock.metier <- names(fleet@metiers[[k]]@catches)
    
    for(k in 1:length(nms.metiers)){
      nms.stock.metier <- names(fleet@metiers[[k]]@catches)
      
      for(j in 1:length(nms.stock.metier)){
       
        #LANDINGS.N,LANDINGS.WT,DISCARDS.N,DISCARDS.WT
        landings.n.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@landings.n)
        landings.n.df$age <- factor(landings.n.df$age)
        landings.n.df$indicator <- 'landings.n'

        discards.n.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@discards.n)
        discards.n.df$age <- factor(discards.n.df$age)
        discards.n.df$indicator <- 'discards.n'

        landings.wt.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@landings.wt)
        landings.wt.df$age <- factor(landings.wt.df$age)
        landings.wt.df$indicator <-'landings.wt'

        discards.wt.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@discards.wt)
        discards.wt.df$age <- factor(discards.wt.df$age)
        discards.wt.df$indicator <- 'discards.wt'
  
        df <- rbind(landings.n.df,discards.n.df,landings.wt.df,discards.wt.df)
        df$stock <- rep(paste(nms.metiers[[k]],"//",nms.stock.metier[j]),dim(df)[1])
        temp <- aggregate(data ~ age + year+stock+indicator, data = df, mean)
        
        p <- ggplot(data=temp, aes(x=year, y=data, fill=age))  + geom_line() +
              geom_point(size=2, shape=21)+
              facet_grid(indicator~stock,scales=c("free_y"))
        print(p)

        alpha.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@alpha)
        alpha.df$age <- factor(alpha.df$age)
        alpha.df$indicator <- 'alpha'

        beta.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@beta)
        beta.df$age <- factor(beta.df$age)
        beta.df$indicator <- 'beta'
        
        catch.q.df <- as.data.frame(fleet@metiers[[k]]@catches[[j]]@catch.q)
        catch.q.df$age <- factor(catch.q.df$age)
        catch.q.df$indicator <- 'catch.q'

        df <- rbind(alpha.df,beta.df,catch.q.df)
        df$stock <- rep(paste(nms.metiers[[k]],"//",nms.stock.metier[j]),dim(df)[1])
        temp <- aggregate(data ~ age + year+stock+indicator, data = df, mean)
        
        p <- ggplot(data=temp, aes(x=year, y=data, fill=age))  + geom_line() +
          geom_point(size=2, shape=21)+
          facet_grid(indicator~stock,scales=c("free_y"))
        print(p)

      }
    }
    dev.off()
  }  
}


