#library(ks)
library(gtools)
library(dismo)

options(warn=-1)

# Retunnnth AUX given the pa and the prob
auc<-function(p0,pa) {
  mean(sample(p0[which(pa==1)],1000,replace=T)>sample(p0[which(pa!=1)],1000,replace=T),na.rm=T)
}

range01 <- function(x) x/sum(x,na.rm=T)
	
#' membership
#' 
#' membership is used to develop an environmental model from known locations, 
#' and use them to predict further occurrences
#' 
#' The main application is species niche modelling.  
#' The algorithm implemented here claims to be both more accurate and 
#' more informative than other methods (hence the name 'WhyWhere').
#'
#' @param files An environmental file
#' @param Pres An array of presence coords 
#' @param Back Absence or background coords
#' @param factor Is a factor variable
#' @param extent The extent of the range of values
#' @return A membership function for that variable


	
membership<-function(x,...) UseMethod("membership")
	
membership<-function(file,P,A,factor=F,extent) {
  ras=raster(file)
  if (strsplit(file,"[.]")[2]=="pgm") {
    xmin(ras)=-180;xmax(ras)=180;ymin(ras)=-90;ymax(ras)=90}
  data=rbind(P,A)
  ras1=crop(ras,extent)
  Pres=extract(ras1,P)
  Back=extract(ras1,data)
  if (factor) breaks=order(unique(c(Pres,Back)))
  else {
    h1=hist(Back,plot=F)
    breaks=h1$breaks
  }
  h2=hist(Pres,breaks=h1$breaks,plot=F)
  prior=c(0,range01(h1$counts))
  post=c(0,range01(h2$counts))
  response=range01(post/prior)
  response[is.nan(response)]<-0
  #model=lm(response~poly(breaks,3))
  model=loess(response~breaks)
  model$membership=data.frame(prior,post,response)
  #plot.membership(list(model))
  #browser()
  pa=c(rep(1,length(Pres)),rep(0,length(Back)))
  values=c(Pres,Back)
  prob=predict(model,newdata=data.frame(breaks=values))
  AUC=auc(prob,pa)
  model$AUC=AUC
  model$name=file
  model$data=data.frame(pa,values,prob)
  model
}

#' expression
#' 
#' expression is used to develop an environmental model from known locations, 
#' and use them to predict further occurrences
#' 
#' The main application is species niche modelling.  
#' The algorithm implemented here claims to be both more accurate and 
#' more informative than other methods (hence the name 'WhyWhere').
#'
#' @param pa A vector of presence/absence
#' @param data A matrix of memberships for each env value 
#' @param dimensions number of terms
#' @return Possibilities arranged in order of probability

expression<-function(pa,x,dimensions=1) {
  cols <- combn( ncol(x), dimensions)
  nams<-apply( cols, 2, 
               function(z) paste(z, collapse='.'))
  
  p0 <- apply( cols, 2, function(z) pmin(x[,z]))
  colnames(p0) <- nams
  p0
  aucs=apply(p0,2,function(x) auc(x,pa))
  aucs=aucs[order(aucs,decreasing=T)] 
}

  
print.membership<-function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
}

#' predict
#' 
#' predict is used to develop an environmental model from known locations, 
#' and use them to predict further occurrences
#' 
#' The main application is species niche modelling.  
#' The algorithm implemented here claims to be both more accurate and 
#' more informative than other methods (hence the name 'WhyWhere').
#'
#' @param pa A vector of presence/absence
#' @param data A matrix of memberships for each env value 
#' @param dimensions number of terms
#' @return Possibilities arranged in order of probability


predict.membership<-function(files,models,extent) {
  m=mapply(function(x,y) {
    ras=raster(x)
    ras0=crop(ras,extent)
    names(ras0)="breaks"
    result=predict(ras0,y)
  },files,models)
  s=stack(m)
  min(s)
}



plot.membership<-function(o) {
  #Needs to handle multiple obs
  side=floor(sqrt(length(o)))
  par(mfcol=c(1,side))
  lapply(o, function(x) {
    plot.ts(x$membership,ylim=c(0,1),xlab="Values",ylab="Response",
            plot.type="single",col=2:4)
    lines(x$fitted,lwd=3,col=4)
    title(x$name)
    })
}
  

simulate.membership<-function(n=100) {
  env=rnorm(n*10)
  data0=data.frame(sample(env,n))
  data=cbind(data0,data0+rnorm(n))
  data=cbind(data,data[,2]+rnorm(n),rnorm(n),rnorm(n))
  pa=rep(0,n)
  pa[which(exp(-abs(data[,1]))>0.25 & exp(-abs(data[,4]))>0.25)]=1
  #pa[which(exp(-abs(data[,4]))>0.75 | exp(-abs(data[,5]))>0.75)]=1
  result=cbind(pa,data)
  #Do OR variables
  colnames(result)=c("pa","vand1","vand1.2","vand1.3","vandor2","vor1")
  result
  }


