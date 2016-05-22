#library(maptools)
library(dismo)
library(data.table)
library(xtable)
library(ggplot2)
 
tfiles=c("s00sd1.25.pgm","fnocwat.pgm")
path1="/home/davids99us/data/Terrestrial"
Tfiles <- list.files(path1,pattern='pgm')
path=paste(system.file(package="dismo"),'/ex',sep='')
files <- list.files(path, pattern='grd')
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
filesb="bio7.grd"
Pres$species=NULL
#browser()
#p=presample(Pres,paste(path,files[1],sep="/"))
result=ww(Pres,Tfiles,dirname=path1,multi=1)
 
plot.ww1(result)
#g=goglm(Pres,result$result$file,dirname=path)
#table1=g[,WW:=result$result$WW]





