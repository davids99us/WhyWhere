#library(maptools)
library(dismo)
library(data.table)


path=paste(system.file(package="dismo"),'/ex',sep='')
files <- list.files(path, pattern='grd')
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
 
o=ww(Pres,files,dirname=path)
 

par(mfcol=c(2,2), mar=c(2, 2, 2, 2) + 0.1)
plot(o$ras,main=paste("Variable",o$name))
ext=c(-86.9333,-39.0667,-24.4500,14.9500)
biome=getraster(paste(path,"biome.grd",sep="/"),ext)
plot(biome)
title("Variable biome")
plot.dseg(o)
l = predict.dseg(o)
plot(l,main="Prediction")
points(o$data$lon,o$data$lat)

#g=goglm(Pres,result$result$file,dirname=path)
#table1=g[,WW:=result$result$WW]
rex1=data.table()
for (i in files) {
e=presample(Pres,paste(path,i,sep="/"))
w0=ww(e,i,dirname=path,type="pa")
auc0=evaluate.ww(e,w0)$auc
auc1=myglm(e)$auc
rex1=rbind(rex1,list(files=i,AUC.WW2=auc0,AUC.GLM=auc1))
}

setorder(rex1,-AUC.WW2)


