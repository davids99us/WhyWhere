#Run on GED file
path="/home/davids99us/data/Terrestrial"
Tfiles=list.files(path,pattern='pgm')
file <- paste(system.file(package="dismo"), '/ex/bradypus.csv',sep='')
Pres <- fread(file,  header=T,sep=",")
Pres$species=NULL
if (FALSE) {
  result=ww(Pres,Tfiles,dirname=path)
save(result,file="resultall6.Rda")
} else load("../resultall6.Rda")

files=result$result$file[1:10]
rex4=data.table()
for (i in files) {
  e=presample(Pres,paste(path,i,sep="/"))
  w0=ww(e,i,dirname=path,type="pa")
  auc0=evaluate.ww(e,w0)$auc
  auc1=myglm(e)$auc
  rex4=rbind(rex4,list(files=i,AUC.WW2=auc0,AUC.GLM=auc1))
}

setorder(rex4,-AUC.WW2)