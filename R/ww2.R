library(data.table)
library(dismo)
library(raster)
library(rgdal)
library(gridExtra)
library(discretization)
library(pROC)
library(ROCR)
library(xtable)
#library(ncdf)
verbose.flag=0

# Get Bradypust test from dismo
file <- paste(system.file(package = "dismo"), "/ex/bradypus.csv", sep = "")
Bradypus_Pres <- fread(file, header = T, sep = ",")
Bradypus_Pres$species = NULL
Bradypus_files <- list.files(path = paste(system.file(package = "dismo"), "/ex", 
    sep = ""), pattern = "grd", full.names = TRUE)
# options(warn=-1)

 
#' range01
#' range01 adjusts range to equal area 
#' @param x vector 
#' @return scaled vector 
range0S <- function(x) x/sum(x, na.rm = T)
range01 <- function(x) (x-min(x))/(max(x)-min(x))

#' getraster
#' getraste gets a raster using raster package, and added 
#' global range for pgm files 
#' @param file file name 
#' @return raster 
#' @return extent c(xmin, xmax, ymin, ymax)
getraster <- function(file,ext) {
    ras = raster(file)
    # In the case is one of my pgm files with no geo information
    if (substr(file, nchar(file) - 3, nchar(file)) == ".pgm") {
        xmin(ras) = -180
        xmax(ras) = 180
        ymin(ras) = -90
        ymax(ras) = 90
        #ras=subs(ras,data.frame(0,NA))
    }
    if (length(ext) == 4) {
    if (ext[1]<xmin(ras) || ext[2]>xmax(ras)) print("Warning X out of range")
    if (ext[3]<ymin(ras) || ext[4]>ymax(ras)) print("Warning Y out of range")
    ras = crop(ras, ext)
    } 
    ras
}

#' getlookup
#' gets the lookup table with probability in each range 
#' @param e the pa values table 
#' @return  lookup table

getlookup<-function(e,ras,cuts) {
#browser()
  lookup = data.table(factors = levels(e$factors), levels = 1:nlevels(e$factors))
setkey(lookup, levels)
# distribution of survey
lookup[, G:=e[pa==0, as.numeric(table(factors))/.N]]
# distribution of cells - dies for large rasters
browser()
if (FALSE) {
  B=table(as.vector(cut(ras,cuts)))
  B=B/sum(B)
  d=data.table(levels=as.numeric(names(B)),B)
  setkey(d,levels)
  lookup=d[lookup]
  lookup[, B:=as.numeric(V2)]
}
#lookup[is.na(get(B)),B:=0,with=FALSE]
#lookup[, B:=e[, as.numeric(table(factors))/.N]]
#distribution of occurrances
lookup[, S:=e[pa == 1, as.numeric(table(factors))/.N]]
lookup
}

#' getlookup
#' gets the lookup table with probability in each range 
#' @param e the pa values table 
#' @return  lookup table

getlookup1<-function(e,ras,type) {
  if(type=="po") {
    lookup=na.omit(as.data.table(freq(ras)))
    temp=e[,.N,by=value]
  }
  if (type=="pa") {
    lookup=e[,.N,by=value]
    setnames(lookup,"N","count")
    temp=e[pa==1,.N,by=value]
  }
  setkey(lookup,value);setkey(temp,value)
  lookup=temp[lookup]
  lookup[is.na(N),N:=0]
  lookup
}


odds<-function(S,G) {
  G[is.na(G)]<-0
  S[is.na(S)]<-0
  S=S/sum(S);G=G/sum(G)
  S=replace(S,S==1,0.99)
  G=replace(G,G==0,0.01)
  (S*(1-G))/(G*(1-S))
}




#' dseg1
#' 
#' dseg1 generates the ww model - a segmented distribution 
#' 
#' @param files A single environmental file
#' @param ext a geographic extent 
#' @param Basis Presence locations and absence locations - pa,  lon, lat - pa a vector of zeros and ones for presence and absence locations
#' @param min.res minimum resolution (xXy)
#' @param min.breaks minimum unique elements
#' @return A model object membership 

dseg1 <- function(file, ext, data, type="po",remove.val=NA,
                  min.res=100,min.breaks=3,
                  max.breaks=10,silent=FALSE,
                  segment="quantile", unique=1) {
  ras = getraster(file,ext)
  data[,value:=extract(ras,cbind(lon,lat))]
  counts=getlookup1(data,ras,type)
   
  #Delete zero row for my pgm files
  if (!is.na(remove.val)) counts=counts[value!=remove.val]
  #browser()
  if (sum(counts$N)==0) return(NULL)
  #Not enought rows
  if (dim(counts)[1]<2) return(NULL)
  #browser()
  #nbreaks = min(max.breaks,dim(counts)[1],floor((length(data$values))^(1/3)))
  nbreaks=10
  freq=1:10/10
  stats=data.table(file=basename(file),XN=dim(ras)[1],
                   YN=dim(ras)[2],breaks=nbreaks)
  stats[,Rbad:=(XN*YN< min.res)][,Bbad:=(nbreaks<min.breaks)][,Ubad:=(dim(counts)[1]<2)]
  verbose(head(stats),"stats")
   
  if (stats$Rbad | stats$Bbad | stats$Ubad ) return(NULL)
  
  bin.size=sum(counts$count)/nbreaks
  r=0;ID=c();j=1;cuts=c(counts$value[1]-1)
  for (i in counts$value) {
    if (r<bin.size) r=r+counts[value==i,count]
    else {r=0;j=j+1;cuts=c(cuts,i)}
    ID=c(ID,j)
  }
  if (length(cuts)<2) return(NULL)
  counts[,ID:=ID]
  width=diff(c(cuts,max(counts$value)))
  lookup <- counts[, lapply(.SD,sum), by=ID]
  lookup[,cuts:=cuts]
  lookup[,width:=width]
  lookup[, odds:=odds(N,count)]
  lookup[, prob:=odds/(1+odds)]
  verbose(lookup,"lookup")
    
  #auc=entropy(lookup$prob)
  if (length(lookup$N)<2 | length(lookup$count)<2) browser()
  #if (file=="/home/davids99us/data/Terrestrial/a00an1.30.pgm") browser()
   
  auc=chisq.test(lookup$N,p=lookup$count,rescale.p=TRUE)$statistic
  o = list(name = names(ras), WW = auc, data = data, raster = ras, 
           lookup = lookup)
  o
}





myglm<-function(e) {
#if (dim(e)[2]==5) {
gm<-glm(pa~value+I(value^2)+I(value^3),data=e, 
        family=binomial(link="logit"))
#gm1<-glm(pa~poly(values,3),data=e, family=binomial(link="logit"))
#} else if (dim(e)[2]==6) {
#  e=na.omit(e)
#gm<-glm(pa~poly(values,3)+poly(V2,3),data=e, family=binomial(link="logit"))
#}
ev<-evaluate(e[pa==1],e,gm)
auc=round(slot(ev,"auc"),3)
list(model=gm,auc=auc,result=ev)
}


#' predict
#' 
#' predict is used to develop an environmental map from known locations, 
#' and use them to predict further occurrences
#' 
#' The main application is species niche modelling.  
#' The algorithm implemented here claims to be both more accurate and 
#' more informative than other methods (hence the name 'WhyWhere').
#'
#' @param obj is a model object from whywhere
#' @return A raster of probability

factor2breaks<-function(x) {
  as.numeric(unique(unlist(strsplit(substr(x, 2, nchar(x) - 1), ","))))
}

predict.dseg <- function(obj) {
  #browser()
     #n1 = factor2breaks(obj$lookup$factors)
     ras = cut(obj$ras, obj$lookup$cuts)
     id = data.frame(cbind(obj$lookup$ID,obj$lookup$prob))
     #browser()
  s = subs(ras, id)
}

#' ww
#' 
#' ww is used to develop an environmental model from known locations, 
#' and use them to predict further occurrences.  
#' The main application is species niche modelling.  
#'
#' @param Pres A 2 column data.frame of known locations
#' @param files The environmental files for raster 
#' @param multi 1 max numer of conjuncts
#' @return result A frame of results
#' @examples
#' result=ww(Bradypus_Pres,Bradypus_files)
#' print(result)



ww <- function(Pres, files, ...) UseMethod("ww")

ww.default <- function(x, y, ...)
  {
    #x <- as.matrix(x)
    #y <- as.numeric(y)
    est <- model.ww(x, y, ...)
    #est$fitted.values <- as.vector(x %*% est$coefficients)
    #est$residuals <- y - est$fitted.values
    #est$call <- match.call()
    class(est) <- "ww"
    est
  }



#' presample
#' presample prepares data for GLM  
#' @param Pres A 2 column data known locations
#' @param files The environmental files for raster  
#' @param trim "crop" around presence points  
#' @param e margin for trim
#' @return data table of locations, pa and values 

presample<-function(data, file=NA,trim=TRUE,e=1,extent=NA) {
  #Get the mask and extent
  if (trim && is.na(extent)) extent=as.numeric(sapply(data[,1:2,with=FALSE], range)) + c(-e, e, -e, e)
  #Get absence points 
  ras=getraster(file[1],extent)
  back=data.table(pa=0,randomPoints(ras,1000))
  setnames(back,c("pa","x","y"),c("pa","lon","lat"))
  pres=data.table(pa=1,lon=data$lon,lat=data$lat)
  data=rbind(pres,back)
  data[,value:=extract(ras,cbind(lon,lat))]
    #data[,eval(paste("V",j,sep="")):=extract(ras,cbind(lon,lat))]
  data[,ID:=.I]
  data
}

goglm<-function(data,files,dirname) {
  Pfiles=paste(path,files,sep="/")
  g=c()
for (i in Pfiles)  {
  #Run GLM for comparison
  e= presample(data,i)
  g=c(g,myglm(e)$auc)
}
data.table(file=files,GLM=g)
}

verbose<-function(data,text,v=NA) {
  if (is.na(v)) v=verbose.flag
  if (v==0) return()
  #browser()
  out=switch(class(data)[1],
  data.table= {
    #extent=as.numeric(sapply(data, range))
    #paste(names(data),extent)
    data
  },
  RasterLayer={
    h=hist(data,breaks=length(unique(data)),plot=FALSE)$density
    if(v<3) length(h) else (paste(h))
  },
  numeric={
    paste(data)
  },
  integer={
    paste(data)
  },
  list={
    paste(list)
  },
  "class not implemented yet")
  print(text)
  print(out)
}

model.ww <- function(data, files, dirname=".", multi=1,plot=FALSE, 
                     limit=0, type="po",remove.val=NA,
                     segment="quantile", trim=TRUE, e=1) {
  result = data.table()
  if (trim) extent=as.numeric(c(range(data$lon),range(data$lat))) + c(-e, e, -e, e)
  else extent=NA
  #print(extent)
  #Determine the type of analysis from data
  #type.code=c("a","po","pa")
  #if (is.na(type)) type=type.code[dim(data)[2]]
   
  #Go through all files
  for (i in files) {
    if (verbose.flag>0) print(i)
    filename=paste(dirname,i,sep="/")
    a = dseg1(filename, extent, data, remove.val=remove.val,type=type)
    if (is.null(a)) {
      print(paste("Skipping",basename(i)))
    } else {
      
    result=rbind(result,list(file=i,WW=a$WW))
    setorder(result, -WW)
    
    #print(paste(a$name,a$AUC))
    if (i == result[1]$file) {
      topdata = a$data #current is best - save pa probs
      if (plot) {
        plot.ww(a)
        #browser()
      }
      #Do the multi variable option if less conjunct than multi
    } else if (multi > length(strsplit(result[1]$file,'&')[[1]])) {
      p = pmin(a$data$prob, topdata$prob)
      #p = a$data$prob * topdata$prob
      auc = auc(data.table(pa=a$data$pa, prob=p))
      r0=list(file=paste(result[1]$file,"&",i,sep =""),WW=auc)
      result=rbind(result,r0)
      setorder(result, -WW)
    }}}
  
  #End of main loop p - get result again
  verbose(head(result,5),"result")
  #browser()
  #print(result)
  if (nrow(result)==0) {
    print("No good files found - returning")
    return(NULL)
  }
  exp=strsplit(result[1]$file,"&")
  filename=paste(dirname,exp[[1]][1],sep="/")
  a = dseg1(filename, extent, data,remove.val=remove.val,type=type)
  a$result = result
  a$dirname=dirname
  #setkey(a$lookup,levels)
  #list(result,dirname,extent,ras) may not need dseg here / do multi
  #browser()
  a
}






#' plot
#' 
#' plot plots the ww model 
#' @param o is a model object from whywhere

plot.ww1<-function(o) {
  par(mfcol=c(2,2), mar=c(2, 2, 2, 2) + 0.1)
  plot(o$ras,main=paste("Variable",o$name))
  plot.dseg(o)
  l = predict.dseg(o)
  plot(l,main="Prediction")
  points(o$data$lon,o$data$lat)
}

plot.dseg<-function(o,legend.loc="topright") {
  mids=barplot(o$lookup$prob,o$lookup$width,
               names.arg=o$lookup$cuts,main="Model")
  legend( legend.loc, legend=c("P(E)","P(E|S)","P(S|E)"),lty=c(1,2,0), 
          pch=c(1,2,15),bg=c(0,0,"gray") )
  prior=o$lookup$count/sum(o$lookup$count)
  post=o$lookup$N/sum(o$lookup$N)
  lines(mids,prior,lty=1)
  lines(mids,post,lty=2)
  points(mids,prior,pch=1)
  points(mids,post,pch=2)
}


#From dismo evaluate
myauc<-function(data) {
  #data=na.omit(data)
  p=data[pa==1,prob]
  a=data[pa==0,prob]
  np=length(p);na=length(a)
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
  auc <- round(R/(na * np),4)
  auc
}


evaluate.ww<-function(test,obj,remove.val=NA) {
  Basis=test[,.(lon,lat)]
  file=obj$result$file[1]
  filename=paste(obj$dirname,obj$result$file[1],sep="/")
  ras = getraster(filename,test$ext)
  e=test[,values := as.numeric(extract(ras, Basis))]  
  
  if (!is.na(remove.val)) e=e[value!=remove.val]
  #browser()
  #Do prediction - cant use predict.dseg as raster based
  cuts=c(obj$lookup$cuts,max(e$value)+1)
  p1=sapply(e$value,function(x) which(cuts>x)[1]-1)
  #browser()
  e[,prob:=obj$lookup$prob[p1]]
  list(auc=myauc(e),fitted.values=e$prob,pa=e$pa)
}



#' dokfold
#' 
#' does the k fold testing
#' @param data_set is a data set from presample
#' @param k is the number of reps
#' @param files is the files

dokfold<-function(data,k=5,files=files,dirname=".",segment="quantile") {
  kfold=1:k
  data_set=as.data.table(data)
  data_set[,k:=kfold]
  rstat=data.table()
  for (i in kfold) {
    train_set=data_set[k!=i]
    train_set<-train_set[,k:=NULL]
    w0=ww(e,files,dirname=dirname,type="pa")
    w0a=evaluate.ww(e,w0)
    g0=glm(pa~values+I(values^2)+I(values^3),data=e, family=binomial(link="logit"))
    ev<-evaluate(e[pa==1],e,g0)
    g0a=round(slot(ev,"auc"),3)
    
    #test set
    test_set=data_set[k==i]
    test_set$k=NULL
    w1a=evaluate.ww(test_set,w0)
    ev=evaluate(test_set[pa==1],test_set,g0)
    g1a=round(slot(ev,"auc"),3)
    
    stat=list(file=w0$result$file[1],GLMtrain=g0a,GLMtest=g1a,
              WWtrain=w0a$auc,WWtest=w1a$auc)
    #browser()
    rstat=rbind(rstat,stat)
  }
  rstat=as.data.frame(na.omit(rstat))
  m=apply(rstat[,2:k],2,mean)
  s=apply(rstat[,2:k],2,sd)
  x=data.table(rbind(m,s/sqrt(k)))
  x[,file:=c("mean","s.e.")]
  rstat=rbind(rstat,x)
}

quanff <- function(vals,freq,quant) {
  ord <- order(vals)
  cs <- cumsum(freq[ord])
  browser()
return(vals[max(which(cs<quant))+1])
}
