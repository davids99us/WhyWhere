library(maptools)
library(dismo)
library(dplyr)

data(wrld_simpl)
files <- list.files(path=paste(system.file(package="dismo"), '/ex',sep=''), 
                    pattern='grd', full.names=TRUE )
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres=as.tbl(read.table(file, header=TRUE,sep=",")) %>% 
  dplyr::select(-matches("species"))
Predictors=stack(files) 

ext=Pres %>% summarise_each(funs(min,max)) + c(-0.5,0.5,-0.5,0.5)
Back=data_frame(lon=runif(100,ext[,1],ext[,3]),lat=runif(100,ext[,2],ext[,4]))
Locations=bind_rows(mutate(Pres,pa=1),mutate(Back,pa=0))
Basis=cbind(Locations,extract(predictors,as.data.frame(Locations[,1:2])))
  
#count(data,biome)
#filter(data,pa==1) %>% count(biome)

p0=Basis %>% group_by(biome) %>%
  summarise(p0=n()) %>% mutate(pp0=p0/sum(p0))
p1=Basis %>% filter(pa==1) %>% group_by(biome) %>%
  summarise(p1=n()) %>% mutate(pp1=p1/sum(p1))
p2 = p0 %>% left_join(p1) %>% 
  mutate(p10=(p1/p0)/sum(p1/p0,na.rm=T)) %>%
  dplyr::select(biome,pp0,pp1,p10)
g=ggplot() + 
  geom_line(aes(biome, p10, colour="blue"), p2) +  
  geom_line(aes(biome, pp1, colour="green"), p1) +
  geom_line(aes(biome, pp0, colour="red"), p0)
print(g)

browser()

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

