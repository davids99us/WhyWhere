library(maptools)
library(dismo)
library(data.table)

#Get presence points
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- read.table(file,  header=T,sep=",")
Pres$species=NULL
#Get predictors
files <- list.files(path=paste(system.file(package="dismo"),'/ex',sep=''), 
                    pattern='grd', full.names=TRUE )

result=whywhere(Pres,files)
print(result)
