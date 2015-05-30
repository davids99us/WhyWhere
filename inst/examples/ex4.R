Tfiles <- list.files(path="/home/davids99us/data/Terrestrial", 
                     pattern='pgm', full.names=TRUE )
train_set=presample(Pres,files[83])
result=ww(train_set,Tfiles,plot=FALSE)
save(result,file="result.Rda")
plot.ww(result)

