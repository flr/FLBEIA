#-------------------------------------------------------------------------------
#           OBSERVATION MODEL FUNCTIONS
#   - fwdBD(stock,fwdControl)  
#
# Dorleta GarcYYYa
# Created: 16/12/2010 07:57:24
# Changed: 16/12/2010 07:57:19
#-------------------------------------------------------------------------------
# fwd.R
# FLash/R/fwd.R
# Copyright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, Cefas
# Last Change: 06 Mar 2009 19:17
# $Id: fwdControl.R 366 2009-10-26 09:37:13Z lauriekell $

trgtNms    <-function() return(c("year","min","val","max","quantity", "rel.year"))
# effNms     <-function() return(c("year","min","val","max","fleet","metier","rel.year","rel.fleet","rel.metier","rel.bound"))
quantityNms<-function() return(c("biomass","catch","landings","discards","f", "f.landings","f.discards"))

validFwdBDControl <- function(object){
	return(TRUE)

  if (dim(object@target)[1]!=dim(object@trgtArray)[1]){
     warning("rows in target & trgtArray don't match")
     return(FALSE)}

  if (any(object@target[,"quantity"] %in% names(quantity))){
     warning("quantity not recognised")
     return(FALSE)}

#  if (length(slot(object, 'effort'))>0){
#     if (dim(object@effort)[1]!=dim(object@effArray)[1]){
#        warning("rows in effort & effArray don't match")
#        return(FALSE)}

#     if (dim(object@target)[1]!=dim(object@effort)[1]){
#        warning("rows in target & effort don't match")
#        return(FALSE)}

#    if (dim(object@trgtArray)[3]!=dim(object@effArray)[3]){
#        warning("iter in trgtArray & effArray don't match")
 #       return(FALSE)}
#     }

	# Everything is fine
	return(TRUE)
  }

setClass("fwdBDControl",
	representation(
		target   ="data.frame",
#		effort   ="data.frame",
        trgtArray="array",
#		effArray ="array",
		block    ="numeric"), ## specifies if mulitple rows are done together
	prototype=prototype(
		target   =data.frame(NULL),
#		effort   =data.frame(NULL),
    trgtArray=array(),
#		effArray =array(),
		block    =numeric()),
	validity=validFwdBDControl
  )

#if (!isGeneric("fwdBDControl")) {
#	setGeneric("fwdBDControl", function(object, ...){
#		value  <-  standardGeneric("fwdBDControl")
#		value
#	})}

#setMethod("fwdBDControl", signature(object="data.frame"),
fwdBDControl <-function(object,trgtArray=NULL,...){

    ##### Internal Functions ###################################################
    setArray<-function(x,nrws,nits=NULL,type="trgtArray"){
       if (is(x,"list") & any(names(x) %in% c("min","val","max"))){
         if (!all(lapply(x,class) %in% c("array","matrix","numeric")))
            stop(paste(type,": elements of list neither 'array', 'matrix' or 'numeric'"))

         if (is.null(nits))
            if      (is(x[[1]],"numeric"))                     nits<-length(x[[1]])
            else if (is(x[[1]],"array") | is(x[[1]],"matrix")) nits<-dim(x[[1]])[length(dim(x[[1]]))]
            else stop("")

         res<-array(NA,dim=c(nrws,3,nits),dimnames=list(1:nrws,c("min","val","max"),iter=1:nits))
         if ("val" %in% names(x)){
            if (is.vector(x$val)) x$val<-array(x$val,dim=c(1,length(x$val)))
            if (nits == dim(x$val)[2])
               res[,"val",]<-x$val
            }
         if ("min" %in% names(x)){
            if (is.vector(x$min)) x$min<-array(x$min,dim=c(1,length(x$min)))
            if (nits == dim(x$min)[2])
               res[,"min",]<-x$min
            }
         if ("max" %in% names(x)) {
            if (is.vector(x$max)) x$max<-array(x$max,dim=c(1,length(x$max)))
            if (nits == dim(x$max)[2])
               res[,"max",]<-x$max}
            }
       else if (is(x,"array") & (length(dim(x))==3)){
          if (is.null(nits))
             nits<-dim(x)[3]

          res<-array(NA,dim=c(nrws,3,nits),dimnames=list(1:nrws,c("min","val","max"),iter=1:nits))

          res[dimnames(x)[[1]],dimnames(x)[[2]],]<-x
          }
       else stop("Has to be either a 3D array or list with 'min', 'max' or 'val' vectors")

       return(res)
       }
	# Creates data.frame with desired column names (nms) and no. of rows (no. yrs)
	# Used for creating target and effort dataframes
    df<-function(yrs,nms){
      df<-NULL
      for (i in nms)
         df<-cbind(df,rep(NA,length(yrs)))
      dimnames(df)<-list(1:length(yrs),nms)
      return(data.frame(df))
      }

    checkMinMax<-function(object)
        {
        # check that if max or min specified then no target & vice versa
        if (any((!is.na(object[,"min"]) | !is.na(object[,"max"])) & !is.na(object[,"val"]))) {
           cat("Can't specify val and both a min or max values")
           return(FALSE)}
        else if (any((!is.na(object[,"min"]) & !is.na(object[,"max"])) & object[,"max"]<=object[,"min"])){
           cat("max less than than min value")
           return(FALSE)}
        else
           return(TRUE)
        }
    ##### End Internal Functions ###############################################

    if (!is(object,"data.frame"))
       stop("target not data.frame")

    if (!("year" %in% names(object)))
       stop("year not specified in object")
    yrs<-object[,"year"]

    res<-new("fwdBDControl")

    ##Targets ##################################################################
    ## Create complete target data frame
    res@target<-df(yrs,trgtNms())
    res@target[,dimnames(object)[[2]]]<-object[,dimnames(object)[[2]]]
    if (!checkBDTarget(res@target))
       stop("target not valid")

    if (!is.null(trgtArray)){
       res@trgtArray<-setArray(trgtArray,length(yrs),type="trgtArray")
       if (length(dim(res@trgtArray[,1,]))==2){
          res@target[,"min"]<-apply(res@trgtArray[,"min",],1,median)
          res@target[,"max"]<-apply(res@trgtArray[,"max",],1,median)
          res@target[,"val"]<-apply(res@trgtArray[,"val",],1,median)}
      else{
          res@target[,"min"]<-median(res@trgtArray[,"min",])
          res@target[,"max"]<-median(res@trgtArray[,"max",])
          res@target[,"val"]<-median(res@trgtArray[,"val",])}}
    else{
       res@trgtArray<-array(as.numeric(NA),dim=c(length(res@target[,1]),3,1),dimnames=list(1:length(res@target[,1]),c("min","val","max"),iter=1))}

    res@target[,"quantity"]<-factor(res@target[,"quantity"],levels=c("ssb","biomass","catch","landings","discards","f","f.landings","f.discards"))

    for (i in 1:length(res@target[,1])){
       if (any(is.na(res@trgtArray[i,"min",]))) res@trgtArray[i,"min",]<-res@target[i,"min"]
       if (any(is.na(res@trgtArray[i,"val",]))) res@trgtArray[i,"val",]<-res@target[i,"val"]
       if (any(is.na(res@trgtArray[i,"max",]))) res@trgtArray[i,"max",]<-res@target[i,"max"]}

    if (!checkMinMax(res@target)) {
       cat(" in target\n")
       stop()}

   
   return(res)
   } # )

showArray<-function(object){
    if(dim(object)[3] > 1){
		  v1 <- apply(object, 1:2, median, na.rm=TRUE)
  		v2 <- apply(object, 1:2, mad,    na.rm=TRUE)
      v3 <- paste(format(v1,digits=5),"(", format(v2, digits=3), ")", sep="")}
    else
      v3 <- paste(format(apply(object, 1:2, median, na.rm=TRUE),digits=5))

    print(array(v3, dim=dim(object)[1:2], dimnames=dimnames(object)[1:2]), quote=FALSE)

		if(dim(object)[3] != 1)
			cat("iter: ", dim(object)[3],"\n\n")}

setMethod('show', signature(object='fwdBDControl'),
  function(object){

  showDFTarget<-function(object){

      nm      <-names(object@target)
      optional<-c("rel.year")
      flag    <-apply(as.matrix(!is.na(object@target[,optional])),2,any)
      
      print(object@target[,c("year","quantity","min","val","max",names(flag[flag]))])

      cat("\n")}



  cat("\nTarget\n")
  showDFTarget(object)
  if (any(!is.na(object@trgtArray)))
     showArray(object@trgtArray)

  })

chkFwdBDControl<-function(ctrl,sr,x,y=NULL){
   if (is(x,"FLStock")){

      return(ctrl)
      }

   }

checkBDTarget<-function(target)
    {
    # check that if max or min specified then no target & vice versa
    if (any((!is.na(target[,"min"]) | !is.na(target[,"max"])) & !is.na(target[,"val"]))) {
       warning("Can't specify a val and a min or max values")
       return(FALSE)}

    if (any((!is.na(target[,"min"]) & !is.na(target[,"max"])) & target[,"max"]<=target[,"min"])){
       warning("max less than than min value")
       return(FALSE)}

	# Should also check quantity

    return(TRUE)
    }

matrixBDTarget <- function(target)
    {
    #reorder columns for C code (???)
    target <- target[,trgtNms()]
    for(i in names(target))
        target[,i] <- as.double(target[,i])

    return(matrix(unlist(target),dim(target)))
    }


chkBDTrgtArrayIters <- function(object,trgtArray,sr)
    {
    if (is(object,'FLlst')) object <- object[[1]]
    # get iterations from trgtArray, stock, SR parameters and SR residuals
    its<-sort(unique(c(length(dimnames(trgtArray)$iter), dims(object)$iter, length(dimnames(sr$params[[1]])$iter), length(dimnames(sr$residuals[[1]])$iter))))
    if (length(its)>2 | (length(its)>1 & its[1]!=1)) stop("iter not 1 or n")
    if (length(its)==2 & length(dimnames(trgtArray)$iter == 1)){
        dmns<-dimnames(trgtArray)
        dmns$iter<-1:its[2]
        trgtArray<-array(trgtArray,dim=unlist(lapply(dmns,length)),dimnames=dmns)}

    return(trgtArray)
    }

# check target quantity is factor and that it is currently implemented
chkBDTargetQuantity <- function(target,object)
    {
    ordDmn<-function(dmn,val){
      tmp       <-1:length(dmn)
      names(tmp)<-dmn

      return(tmp[ac(val)])
      }

    if (!is(target[,"quantity"],"factor"))
        target[,"quantity"]<-factor(target[,"quantity"],quantityNms())
    if (!all(as.character(target[,"quantity"]) %in% quantityNms()))
        stop("invalid quantity in control target")

	  return(target)
    }


    