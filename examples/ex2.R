#Get presence points
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
print(file)
Pres <- read.table(file,  header=T,sep=",")
Pres$species=NULL
#Get predictors
files <- list.files(path=paste(system.file(package="dismo"),'/ex',sep=''), 
                    pattern='grd', full.names=TRUE )
Tfiles <- list.files(path="/home/davids99us/Dropbox/R/Temp/WhyWhere/Terrestrial", 
                     pattern='pgm', full.names=TRUE )
#files=files[1:(length(files)-1)] #Remove the factor variable - not implemented
#files="/usr/lib64/R/library/dismo/ex/biome.grd"
#files=c(files,Tfiles)
#files="/home/davids99us/Dropbox/R/Temp/WhyWhere/Terrestrial/a00sd1.33.pgm"

result=whywhere(Pres,files)
print(result)