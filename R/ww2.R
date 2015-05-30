library(data.table)
library(dismo)
library(raster)
library(rgdal)
library(ggplot2)
library(gridExtra)
library(discretization)

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

getlookup<-function(e) {
lookup = data.table(factors = levels(e$factors), levels = 1:nlevels(e$factors))
setkey(lookup, levels)
lookup[, g1:=e[pa==0, table(factors)/.N]]
lookup[, g1:=e[, table(factors)/.N]]
lookup[, g2:=e[pa == 1, table(factors)/.N]]
lookup[, odds:=as.numeric(g2/g1)]
#lookup$odds[is.nan(lookup$odds)] <- 0
lookup[, prob:=odds/(1+odds)]
#lookup[, prob:=log(1/(1+exp(-odds)))]
#lookup[,prob:=range0S(cond*g2)]
#print(lookup)
#browser()
lookup
}

myGLM<-function(e) {
  gm1<-glm(pa~values+I(values^2)+I(values^3),data=e, family=binomial(link="logit"))
  ev1<-evaluate(e[pa==1],e,gm1)
  round(slot(ev1,"auc"),3)
}


#' dseg
#' 
#' dseg generates the ww model - a segmented distribution 
#' 
#' @param files A single environmental file
#' @param ext a geographic extent 
#' @param Basis Presence locations and absence locations - pa,  lon, lat - pa a vector of zeros and ones for presence and absence locations
#' @param min.res minimum resolution (xXy)
#' @param min.breaks minimum unique elements
#' @return A model object membership 

dseg <- function(file, ext, Basis, categorical=FALSE, 
                 plot=FALSE,min.res=100,min.breaks=3, 
                 type=4, segment="quantile", unique=1) {
  pa=Basis$pa
  Basis=Basis[,.(lon,lat)]
  ras = getraster(file,ext)
  if (substr(basename(file),1,4)=="fact") categorical=TRUE
  # Table of data - what about NA's
  e = data.table(pa, Basis, values = extract(ras, Basis))
  #print(paste(length(unique(e$values)),file))
  #browser()
  if (length(unique(e[pa==1]$values))<unique | length(unique(e[pa==0]$values))<unique ) {
    print(paste("Insufficient unique points in",file))
    return(NULL)
  }
  #Run GLM for comparison
  gm1<-glm(pa~values+I(values^2)+I(values^3),data=e, family=binomial(link="logit"))
  #gm1<-glm(pa~poly(values,3),data=e, family=binomial(link="logit"))
  ev1<-evaluate(e[pa==1],e,gm1)
  BAUC=round(slot(ev1,"auc"),3)
  
  #browser()
  cuts=switch(segment,
    entropy={
     r=range(e$values)
     c(r[1],cutPoints(e$values,e$pa),r[2])
    },
    quantile={
      # Values into factor levels using Freedman-Diaconis rule
      nbreaks = min(length(unique(e$values)), floor((length(e$values))^(1/3)))
      stats=data.table(name=basename(file),XN=dim(ras)[1],YN=dim(ras)[2],breaks=nbreaks)
      stats[,N:=(XN*YN< min.res)][,B:=(nbreaks<min.breaks)]
      if (stats$N | stats$B) return(NULL)
        # else print(stats)
      quants = quantile(na.omit(e[pa==0,values]), seq(0, 1, 1/nbreaks), type=type)
    },
    even={
      cuts = min(length(unique(e$values)), floor((length(e$values))^(1/3)))
      #browser()
    },
    categorical=as.factor(e$values))
    #browser()
    factors = try(cut(e$values, cuts), silent=TRUE)
    # If quantiles dont work set values as factors#
    # Usually because too highly skewd - 
    # then breaks are not unique
    if (class(factors) == "try-error") {
    #e[, factors:=as.factor(e$values)]
    e[, factors:=cut(e$values,c(min(e$values)-1,unique(e$values)))]
    categorical=TRUE
    } else e[, factors:=factors]
  setkey(e, factors)
  # Construct lookup table for probs on factors
  lookup=getlookup(e)
 
  #setorder(lookup, factors)
  e = lookup[e]
  auc = auc(e)
   
  o = list(name = names(ras), WW = auc, data = e, raster = ras, 
           lookup = lookup, categorical=categorical, GLM=BAUC, glm=gm1)
  if (plot) {
    plot.ww(o)
    browser()
  }
  o
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
     n1 = factor2breaks(obj$lookup$factors)
     ras = cut(obj$ras, n1)
     id = data.frame(cbind(1:(nrow(obj$lookup)),obj$lookup$prob))
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



#' model.ww
#' model.ww develops the model 
#' @param Pres A 2 column data.frame of known locations
#' @param files The environmental files for raster  
#' @param multi multiple dimensions
#' @param limit minimum for entry into results
#' @param beam number of entries to keep in table
#' @param extent range of plot.  May be vector length 4 to specify.
#' @param trim "crop" around presence points  
#' @param plot logical 
#' @param e margin for trim
#' @return scaled vector 

presample<-function(Pres, mask=NA,trim=TRUE,e=1,split=FALSE) {
  #Get the mask and extent
  if (is.na(mask)) mask=files[1]
  ras=getraster(mask,NA)
  if (length(extent) == 4) ext=extent
  else ext=c(xmin(ras),xmax(ras),ymin(ras),ymax(ras))
  if (trim) ext=as.numeric(sapply(Pres, range)) + c(-e, e, -e, e)
  #Get data
  ras=crop(ras,ext)
  back=data.table(pa=0,randomPoints(ras,1000))
  colnames(back)=c("pa","lon","lat")
  #if (split) pres=cbind(pa=1,Pres[sample(.N,.N/2,replace=TRUE)])
  #else pres=cbind(pa=1,Pres[sample(.N,1000,replace=TRUE)])
  pres=cbind(pa=1,Pres)
  #Combine
  data=rbind(pres,back)
  list(data=data,extent=ext)
}


model.ww <- function(set, files, multi=FALSE,plot=FALSE, limit=0, segment="quantile") {
  result = data.table(name = "limit", WW = limit, GLM=limit,file = NA)
  #print(segment)
  #browser()
  #Go through all files
  for (i in files) {
    a = dseg(i, set$extent, set$data, plot=plot, segment=segment)
    if (is.null(a)) {
      print(paste("Skipping",basename(i)))
    } else {
    result=rbind(result,list(name = a$name, WW=a$WW, GLM=a$GLM, file = i))
    setorder(result, -WW)
    #print(paste(a$name,a$AUC))
    if (a$name == result[1]$name) {
      topdata = a$data #current is best - save pa probs
      if (plot) {
        plot.ww(a)
        #browser()
      }
      #Do the multi variable option
    } else if (multi) {
      p = pmin(a$data$prob, topdata$prob)
      #print(paste(length(a$data$prob), length(topdata$prob), length(p)))
      auc = auc(data.table(pa=a$data$pa, prob=p))
      #browser()
      #Do the glm
      y=na.omit(data.table(pa=a$data$pa, y1=topdata$values, y2=a$data$values))
      #print(file) 
      #print(nrow(y))
      gm2<-glm(pa~poly(y1,3)+poly(y2,3),data=y, family=binomial(link="logit"))
      ev2<-evaluate(y[pa==1],y,gm2)
      glmauc=round(slot(ev2,"auc"),3)
      result=rbind(result,
                 list(name=paste(result[1]$name,".", a$name, sep = ""),
                      WW=auc,GLM=glmauc,file = i))
      setorder(result, -WW)
    }}}
  
  #print(result)
  if (nrow(result)==1) {
    print("No good files found - returning")
    return(NULL)
  }
  a = dseg(result[1]$file, set$ext, set$data)
  a$result = result
  setkey(a$lookup,levels)
  a
}




#' plot
#' 
#' plot plots the ww model 
#' @param o is a model object from whywhere

plot.ww <- function(o,legend="none") {
  a=b=c=d=0
  p0=plot.map(o$ras)+ theme(legend.position=legend,
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(), 
                            axis.text.y = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.ticks.y = element_blank(),
                            plot.margin = unit(c(a,b,c, d),"cm"))
  
  p1=plot.dseg(o) + theme(legend.position=legend)
  #                       axis.text.x=element_blank(),
  #                       axis.ticks.x=element_blank(),
  #                       axis.ticks.y = element_blank(),
  #                       plot.margin = unit(c(a,b,c, d),"cm"))
  l = predict.dseg(o)
  Pres = o$data[pa == 1, .(lon, lat)]
  
  p2=plot.map(l,Pres) + theme(legend.position=legend,
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.text.x = element_blank(), 
                              axis.text.y = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.ticks.y = element_blank(),
                              plot.margin = unit(c(a,b,c, d),"cm"))
  
  
  text2=data.table(Value=c(o$name,o$WW,o$GLM))
  rownames(text2)=c("Variable","WW","GLM")
  
  blankmap <- ggplot(text2,aes(1,1))+geom_blank(aes(3,3))+
    theme(
      plot.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()) +
    theme(plot.margin = unit(c(a,b,c, d),"cm")) +
    annotation_custom(grob = tableGrob(text2),xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  #browser()
  grid.arrange(p0, blankmap, p1, p2, ncol=2) 
}

#' plot.dseg
#' 
#' plot.dseg plots the segregated model
#' @param o is a model object from whywhere

plot.dseg <- function(o) {
  g1=as.numeric(o$lookup$g1)
  g2=as.numeric(o$lookup$g2)
  g1 = data.frame(x=o$lookup$levels,y=g1/sum(g1))
  g2 = data.frame(x=o$lookup$levels,y=g2/sum(g2))
  #print(o$lookup)
  #print(o$name)
  bins=factor2breaks(o$lookup$factors)
  o$lookup[,width:=diff(bins)]
  o$lookup[,Values:=bins[-1]-width/2]

p=ggplot(o$lookup, aes(x=Values, y=prob, width=width)) + 
  geom_bar(aes(fill=1), stat="identity", position="identity") 
  #geom_line(data=o$lookup, aes(x=Values, y=g1), colour="blue") +
  #geom_line(data=o$lookup, aes(x=Values, y=g2), colour="green")
  p
}

plot.map<-function(l,Pres=NULL) {
  #browser()
  map.p <- rasterToPoints(l)
  df <- data.frame(map.p)
  #Make appropriate column headings
  colnames(df) <- c("Longitude", "Latitude", "MAP")
  #Call in point data, in this case a fake transect (csv file with lat and lon coordinates)
  #sites <- data.frame(read.csv(“/your/path/to/pointfile.csv”))
  sites=Pres
  #Now make the map
  
  p2=ggplot(data=df, aes(y=Latitude, x=Longitude)) +
    geom_raster(aes(fill=MAP)) +
    coord_equal(ratio = 1)
  
  if (!is.null(Pres)) p2=p2 + geom_point(data=sites, aes(x=lon, y=lat), color="green", size=3, shape=1) +
    #theme_bw() +
    #coord_equal() +
    #scale_fill_gradient("MAP (prob)", limits=c(0,2500)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16, angle=90),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
    )
  p2
}

#From dismo evaluate
auc<-function(data) {
  data=na.omit(data)
  p=data[pa==1,prob]
  a=data[pa==0,prob]
  np=length(p);na=length(a)
  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
  auc <- round(R/(na * np),4)
  auc
}


evaluate.ww<-function(test,obj) {
  Basis=test$data[,.(lon,lat)]
  ras = getraster(obj$result$file[1],test$ext)
  
  e=test$data[,values := as.numeric(extract(ras, Basis))]  
  #Do the glm  
  ev1<-evaluate(e[pa==1],e[pa==0],obj$glm)
  BAUC=round(slot(ev1,"auc"),2)
  
  #Do prediction - cant use predict.dseg as raster based
  n1 = factor2breaks(obj$lookup$factors)
  e[,factors:=cut(values, n1)]
  setkey(e, factors)
 
    p = obj$lookup[e]
    auc = auc(p)
    o = list(WW = auc, GLM=BAUC)
    o
}


#' dokfold
#' 
#' does the k fold testing
#' @param data_set is a data set from presample
#' @param k is the number of reps
#' @param files is the files

dokfold<-function(data_set,k=5,files=files,segment="quantile") {
  kfold=1:k
  data_set$data[,k:=kfold]
  rstat=data.table(name=NA,WWtr=NA,GLMtr=NA,WWte=NA,GLMte=NA)
  for (i in kfold) {
    train_set=list(data=data_set$data[k!=i],ext=data_set$extent)
    r0=ww(train_set,files,segment=segment)
    test_set=list(data=data_set$data[k==i],ext=data_set$extent)
    r1=evaluate.ww(test_set,r0)
    #browser()
    stat=list(name=r0$name,WWtr=r0$WW,GLMtr=r0$GLM,WWte=r1$WW,GLMte=r1$GLM)
    rstat=rbind(rstat,stat)
  }
  rstat=as.data.frame(na.omit(rstat))
  m=apply(rstat[,2:k],2,mean)
  s=apply(rstat[,2:k],2,sd)
  x=data.table(rbind(m,s))
  x[,name:=c("mean","sd")]
  rstat=rbind(rstat,x)
}

