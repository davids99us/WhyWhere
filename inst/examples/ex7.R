#library(maptools)
library(dismo)
library(data.table)
library(xtable)
library(ggplot2)

par(mfrow=c(2,2))
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL

plot2<-function(r1,r2,title) {
  p1 = prediction(r1$fitted.values,r1$data$pa) 
  p2 = prediction(r2$fitted.values,r2$pa)
  perf1 <- performance(p1, "tpr", "fpr")
  perf2 <- performance(p2, "tpr", "fpr")
  plot(perf1,main=title,col=2)
  lines(perf2@x.values[[1]], perf2@y.values[[1]], col = 1)
  c(performance(p1,'auc'),performance(p2,"auc"))
}
 

path=paste(system.file(package="dismo"),'/ex',sep='')
result=ww(Pres,c("bio7.grd"),dirname=path,plot=FALSE,multi=FALSE)
e=presample(Pres,paste(path,"bio7.grd",sep="/"))
we=evaluate.ww(e,result)
r=myglm(e)
plot2(r$model,we,"bio7.grd")

path=paste(system.file(package="dismo"),'/ex',sep='')
result=ww(Pres,c("biome.grd"),dirname=path,plot=FALSE,multi=FALSE)
e=presample(Pres,paste(path,"biome.grd",sep="/"))
we=evaluate.ww(e,result)
r=myglm(e)
plot2(r$model,we,"biome.grd")


path="/home/davids99us/data/Terrestrial"
result=ww(Pres,c("lccld08.pgm"),dirname=path,plot=FALSE,multi=FALSE)
e=presample(Pres,paste(path,"lccld08.pgm",sep="/"))
we=evaluate.ww(e,result)
r=myglm(e)
plot2(r$model,we,"lccld08.pgm")

path="/home/davids99us/data/Terrestrial"
result=ww(Pres,c("fnocwat.pgm"),dirname=path,plot=FALSE,multi=FALSE)
e=presample(Pres,paste(path,"fnocwat.pgm",sep="/"))
we=evaluate.ww(e,result)
r=myglm(e)
plot2(r$model,we,"fnocwat.pgm")


#path=paste(system.file(package="dismo"),'/ex',sep='')
#result=ww(Pres,c("bio7.grd","biome.grd"),dirname=path,multi=TRUE)
#e=presample(Pres,paste(path,c("bio7.grd","biome.grd"),sep="/"))
#r=myglm(e)
#p=plot2(r$model,result,"bio7.grd&biome.grd")
#path="/home/davids99us/data/Terrestrial"
#result=ww(Pres,c("fnocwat.pgm","i00an1.1.pgm"),dirname=path,plot=FALSE,multi=FALSE)
#e=presample(Pres,paste(path,c("fnocwat.pgm","i00an1.1.pgm"),sep="/"))
#r=myglm(e)
#plot2(r$model,result,"fnocwat.pgm&i00an1.1.pgm")

