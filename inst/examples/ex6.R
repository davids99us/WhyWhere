library(maptools)
library(dismo)
library(data.table)
library(xtable)
#pdf("lerista3.pdf")
path1="/home/davids99us/data/Terrestrial"
Tfiles <- list.files(path1,pattern='pgm')
path2="../data/Lerista"
files <- list.files(path2,pattern='GTiff')
afiles<-c("mgv8610.pgm","mgv8511.pgm","mgv8806.pgm")  
file <- '../data/Lerista/locations.txt'
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file)
Pres$no=NULL
Pres$species=NULL

result=ww(Pres,Tfiles,dirname=path1)
l = predict.dseg(result)

plot.ww1(result)

 