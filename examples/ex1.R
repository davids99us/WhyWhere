library(maptools)
library(dismo)
data(wrld_simpl)
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
bradypus <- read.table(file, header=TRUE,sep=",")
files <- list.files(path=paste(system.file(package="dismo"), '/ex',sep=''), pattern='grd', full.names=TRUE )
Pres=bradypus[,-1]
ext1=c(range(Pres[,1]),range(Pres[,2]))
ext=ext1+c(-0.5,0.5,-0.5,0.5)
Back=cbind(lon=runif(100,ext[1],ext[2]),lat=runif(100,ext[3],ext[4]))
           
models=lapply(files,function(x) membership(x,Pres,Back,extent=ext))

data=sapply(models,function(x) x$data$prob)
cause=expression(models[[1]]$data$pa,data,dimensions=2)

nos=as.numeric(strsplit(names(cause)[1],"[.]")[[1]])
newfiles=files[nos]
Bmodels=models[nos]
result=predict.membership(newfiles,Bmodels,extent=ext)
par(mfcol=c(1,1))
plot(result)
points(Pres, col='blue')

