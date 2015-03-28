library(maptools)
library(dismo)
library(data.table)

#Get presence points
#data(wrld_simpl)
files <- list.files(path=paste(system.file(package="dismo"), 
                               '/ex',sep=''), pattern='grd', full.names=TRUE )
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- read.table(file,  header=T,sep=",")
Pres$species=NULL
result=whywhere(Pres,multi=F,files)
plot.model(result)
plot.membership(result)
