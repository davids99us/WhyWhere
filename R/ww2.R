library(data.table)
library(dismo)
library(raster)
library(rgdal)
library(ggplot2)
library(gridExtra)


# Get Bradypust test from dismo
file <- paste(system.file(package = "dismo"), "/ex/bradypus.csv", sep = "")
Bradypus_Pres <- read.table(file, header = T, sep = ",")
Bradypus_Pres$species = NULL
Bradypus_files <- list.files(path = paste(system.file(package = "dismo"), "/ex", 
    sep = ""), pattern = "grd", full.names = TRUE)
# options(warn=-1)

 
#' range01
#' range01 adjusts range to equal area 
#' @param x vector 
#' @return scaled vector 
range01 <- function(x) x/sum(x, na.rm = T)

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

dseg <- function(file, ext, Basis, categorical=FALSE, plot=FALSE,min.res=100,min.breaks=3) {
  pa=Basis$pa
  Basis=Basis[,.(lon,lat)]
  ras = getraster(file,ext)
  if (substr(basename(file),1,4)=="fact") categorical=TRUE
  # Table of data - what abou NA's
  e = data.table(pa, Basis, values = extract(ras, Basis))
  # Values into factor levels using Freedman-Diaconis rule
  nbreaks = min(length(unique(e$values)), floor((length(e$values))^(1/3)))
  
  #Determine if should apply 
  stats=data.table(name=basename(file),XN=dim(ras)[1],YN=dim(ras)[2],breaks=nbreaks)
  stats[,N:=(XN*YN< min.res)][,B:=(nbreaks<min.breaks)]
  if (stats$N | stats$B) return(NULL)
  # else print(stats)
  
  # Determine factors - if quantiles doesnt work do even breaks
  #browser()
  quants = quantile(na.omit(e$values), seq(0, 1, 1/nbreaks))
  if (categorical) {
    factors=as.factor(e$values)
    #browser()
  } else {
    factors = try(cut(e$values, quants), silent=TRUE)
  }
  # If quantiles dont work set values as factors#
  # Usually because too highly skewd - 
  # then breaks are not unique
  if (class(factors) == "try-error") {
    e[, factors:=as.factor(e$values)]
    categorical=TRUE
  } else e[, factors:=factors]
  setkey(e, factors)
  # Construct lookup table for probs on factors
  lookup = data.table(factors = levels(e$factors), levels = 1:nlevels(e$factors))
  setkey(lookup, factors)
  lookup[, g1:=e[, table(factors)]]
  lookup[, g2:=e[pa == 1, table(factors)]]
  lookup[, prob:=as.numeric(round(g2/g1, 2))]
  lookup$prob[is.nan(lookup$prob)] <- 0
  #setorder(lookup, factors)
  e = lookup[e]
  auc = auc(e)
  o = list(name = names(ras), AUC = auc, data = e, raster = ras, 
           lookup = lookup, categorical=categorical)
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


predict.dseg <- function(obj) {
   if (!obj$categorical) {
     n1 = as.numeric(unique(unlist(
     strsplit(substr(obj$lookup$factors, 2, nchar(obj$lookup$factors) - 1), ","))))
     ras = cut(obj$ras, n1)
     id = data.frame(cbind(1:(nrow(obj$lookup)),obj$lookup$prob))
   }
  else {
    n1 = as.numeric(obj$lookup$factors)
    ras=obj$ras
    id = data.frame(cbind(n1,obj$lookup$prob))
    }
  #m = n1[order(n1)]
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

getSet<-function(Pres,ras,split) {
back=data.table(pa=0,randomPoints(ras,1000))
colnames(back)=c("pa","lon","lat")
if (split) pres=cbind(pa=1,Pres[sample(.N,.N/2,replace=TRUE)])
else pres=cbind(pa=1,Pres[sample(.N,1000,replace=TRUE)])
#Combine
set=rbind(pres,back)
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

model.ww <- function(Pres, files, mask=NA, multi=F,limit=0, 
                     split=FALSE,trim=TRUE,beam=5,extent=NA,plot=FALSE,e=1) {
  #Get the mask and extent
  e=1
  if (is.na(mask)) mask=files[1]
  ras=getraster(mask,NA)
  if (length(extent) == 4) ext=extent
  else ext=c(xmin(ras),xmax(ras),ymin(ras),ymax(ras))
  if (trim) ext=as.numeric(sapply(Pres, range)) + c(-e, e, -e, e)
  #Get data
  ras=crop(ras,ext)
  train_set=getSet(Pres,ras,split)
  if (split) test_set=getSet(Pres,ras,split)
  else test_set = train_set
  
  result = data.table(name = "limit", AUC = limit, file = NA)
  #browser()
  #Go through all files
  for (i in files) {
    a = dseg(i, ext, train_set, plot=FALSE)
    if (is.null(a)) {
      #print(paste("Skipping",basename(i)))
    } else {
    result=rbind(result,list(name = a$name, AUC = a$AUC,file = i))
    setorder(result, -AUC)
    #print(paste(a$name,a$AUC))
    if (a$name == result[1]$name) {
      topdata = a$data #current is best - save pa probs
      if (plot) plot.ww(a)
      #print(result[1:3])
    } else if (multi) {
      p = pmin(a$data$prob, topdata$prob)
      #print(paste(length(a$data$prob), length(topdata$prob), length(p)))
      auc = auc(data.table(pa=a$data$pa, prob=p))
      result=rbind(result,
                 list(name=paste(result[1]$name,".", a$name, sep = ""),
                      AUC=auc,file = i))
      setorder(result, -AUC)
    }}}
  
  #print(result)
  if (nrow(result)==1) {
    print("No good files found - returning")
    return(NULL)
  }
  a = dseg(result[1]$file, ext, train_set)
  a$result = result
  setkey(a$lookup,levels)
  a
}




#' plot
#' 
#' plot plots the ww model 
#' @param o is a model object from whywhere

plot.ww <- function(o) {
  a=b=c=d=0
  p0=plot.map(o$ras)+ theme(legend.position="right",
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(), 
                            axis.text.y = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.ticks.y = element_blank(),
                            plot.margin = unit(c(a,b,c, d),"cm"))
  
  p1=plot.dseg(o) + theme(legend.position="none",
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.ticks.y = element_blank(),
                          plot.margin = unit(c(a,b,c, d),"cm"))
  l = predict.dseg(o)
  Pres = o$data[pa == 1, .(lon, lat)]
  
  p2=plot.map(l,Pres) + theme(legend.position="right",
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.text.x = element_blank(), 
                              axis.text.y = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.ticks.y = element_blank(),
                              plot.margin = unit(c(a,b,c, d),"cm"))
  
  
  text2=data.table(Value=c(o$name,o$AUC))
  rownames(text2)=c("Variable","AUC")
  
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
  
  grid.arrange(p0, blankmap, p1, p2, ncol=2, ncol=2) 
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
  
  p <- ggplot(o$lookup, aes(x=levels, y=prob)) +
    geom_bar(stat="identity", alpha=0.75) +
    geom_line(data=g1, aes(x=x, y=y), colour="blue") +
    geom_line(data=g2, aes(x=x, y=y), colour="green")
  p
}

plot.map<-function(l,Pres=NULL) {
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
  auc <- round(R/(na * np),2)
  auc
}
