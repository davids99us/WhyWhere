#Get presence points
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
print(file)
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
Pres$lat=Pres$lat
#Get predictors
files <- list.files(path=paste(system.file(package="dismo"),'/ex',sep=''), 
                    pattern='grd', full.names=TRUE )
Tfiles <- list.files(path="/home/davids99us/data/Terrestrial", 
                     pattern='pgm', full.names=TRUE )
#files=files[1:(length(files)-1)] #Remove the factor variable - not implemented
#files="/usr/lib64/R/library/dismo/ex/biome.grd"
#files=c(files,Tfiles)
#files="/home/davids99us/Dropbox/R/Temp/WhyWhere/Terrestrial/a00sd1.33.pgm"

result=ww(Pres,Tfiles,plot=TRUE,trim=TRUE,multi=TRUE)
#plot(result)
#print(result)