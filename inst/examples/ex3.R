result1=dokfold(data_set,k=5,files=files,segment="even")
result2=dokfold(data_set,k=5,files=files,segment="quantile")
result3=dokfold(data_set,k=5,files=files,segment="entropy")
r=rbind(result1[6],result2[6],result3[6])
r$file=NULL
r[,segment:=c("even","quantile","entropy")]
