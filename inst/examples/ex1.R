library(maptools)
library(dismo)
library(data.table)
library(xtable)

files <- list.files(path=paste(system.file(package="dismo"), 
                               '/ex',sep=''), pattern='grd', full.names=TRUE )
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
train_set=presample(Pres,files[1])
result=ww(train_set,files,multi=FALSE)
plot.ww(result)



