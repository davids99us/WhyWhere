library(maptools)
library(dismo)
library(data.table)

#Does 5 fold validation
path=paste(system.file(package="dismo"),'/ex',sep='')
files <- list.files(path, pattern='grd')
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
k=5
e=presample(Pres,paste(path,"bio7.grd",sep="/"))
rex2=dokfold(e,k=5,files=c("bio7"),dirname=path)






