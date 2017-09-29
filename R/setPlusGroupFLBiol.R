
#-------------------------------------------------------------------------------
#                           setPlusGroupFLBiol
#   
# Agurtzane Urtizberea
#-------------------------------------------------------------------------------


setPlusGroupFLBiol<-	function(x, plusgroup, na.rm=FALSE)
{
  pg.wt.mean <-c("wt","m","spwn","fec","mat")
  
  #check plusgroup valid
  if (!missing(plusgroup))
    x@range["plusgroup"]<-plusgroup
  if(x@range["plusgroup"] > x@range["max"])
    return("Error : plus group greater than oldest age")
  
  #Perform +grp calcs
  pg.range <- as.character(x@range["max"]:x@range["plusgroup"])
  
  #do the weighted stuff first
  for (i in pg.wt.mean[-c(4,5)]){
    if (dims(n(x))$iter!=dims(slot(x,i))$iter) 
      slot(x,i)<-propagate(slot(x,i),dims(n(x))$iter) 
    slot(x,i)[as.character(x@range["plusgroup"])]<-quantSums(slot(x,i)[pg.range]*x@n[pg.range])/quantSums(x@n[pg.range])
  }
 # i <- "fec"
  if (dims(n(x))$iter!=dims(fec(x))$iter) 
   fec(x)<-propagate(fec(x),dims(n(x))$iter)
  fec(x)[as.character(x@range["plusgroup"])]<-quantSums(fec(x)[pg.range]*x@n[pg.range])/quantSums(x@n[pg.range])
  
  # i <- "mat"
  if (dims(n(x))$iter!=dims(mat(x))$iter) 
    mat(x)<-propagate(mat(x),dims(n(x))$iter)
  mat(x)[as.character(x@range["plusgroup"])]<-quantSums(mat(x)[pg.range]*x@n[pg.range])/quantSums(x@n[pg.range])
  
  x@n[as.character(x@range["plusgroup"])]<-quantSums(x@n[pg.range])
  
  x<-x[as.character(x@range["min"]:x@range["plusgroup"])]
  fec(x) <- fec(x)[as.character(x@range["min"]:x@range["plusgroup"])]
  mat(x) <- mat(x)[as.character(x@range["min"]:x@range["plusgroup"])]
  
  x@range["max"]<-x@range["plusgroup"]
  x@range["minfbar"] <- ifelse( x@range["minfbar"]>x@range["plusgroup"], x@range["plusgroup"], x@range["minfbar"])
  x@range["maxfbar"] <- ifelse( x@range["maxfbar"]>x@range["plusgroup"], x@range["plusgroup"], x@range["maxfbar"])
  
  return(x)
}
