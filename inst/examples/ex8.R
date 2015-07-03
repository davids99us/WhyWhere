#Compare some of the top variables
par(mfrow=c(2,2), mar=c(2, 2, 2, 2) + 0.1)
path="/home/davids99us/data/Terrestrial"
Tfiles=list.files(path,pattern='pgm')
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL

o =ww(Pres,c("lccld08.pgm"),dirname=path)
plot.dseg(o)  
plot(predict.dseg(o))
points(o$data$lon,o$data$lat)
title("lccld07")

o =ww(Pres,c("fnocwat.pgm"),dirname=path)
plot.dseg(o)  
plot(predict.dseg(o))
points(o$data$lon,o$data$lat)
title("fnocwat")
