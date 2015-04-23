library(maptools)
library(dismo)
library(data.table)

#Get presence points
#data(wrld_simpl)
#files <- list.files(path=paste(system.file(package="dismo"), 
#                               '/ex',sep=''), pattern='grd', full.names=TRUE )
files <- list.files(path="/home/davids99us/data/dismo",pattern='grd', full.names=TRUE )
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
result=ww(Pres,files,plot=FALSE,multi=TRUE,beam=100,trim=FALSE)
plot.ww(result)



