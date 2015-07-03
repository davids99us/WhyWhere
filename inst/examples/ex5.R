path=paste(system.file(package="dismo"),'/ex',sep='')
files <- list.files(path, pattern='grd')
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
result=ww(Pres,files,dirname=path,multi=TRUE,plot=FALSE)
plot.ww(result)