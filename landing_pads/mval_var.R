library(data.table)
library(ggplot2)
setwd("/data/borth/hdhiman/janssen/integration_sites/filter_regions/src2/variant_calls/bg")

list<-read.table("/data/borth/hdhiman/janssen/integration_sites/filter_regions/src2/list_var_bg.txt", header=FALSE)
data<-lapply(as.character(paste(list$V1,"bedGraph", sep=".")), fread, header=FALSE)
names(data)<-list$V1
for(i in list$V1)
{
print(i)
ip=paste(i,"bedGraph", sep="."); 
out=paste(i,"mval.bedGraph", sep="_"); 
mad<-median(abs(data[[i]]$V4-median(data[[i]]$V4)))
med<-median(data[[i]]$V4)
print(mad)
print(med)
data[[i]]$mval<-(data[[i]]$V4-median(data[[i]]$V4))/3;
fwrite(subset(data[[i]], select=-c(4)), file=out,sep="\t",quote=F,row.names=F,col.names=F)
}
