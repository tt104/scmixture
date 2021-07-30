library(umap)
library(ggplot2)

# Original data
D<-as.matrix(read.csv(snakemake@input[[1]],header=TRUE,row.names=1))

clusters<-as.matrix(read.csv(snakemake@input[[2]],header=TRUE,row.names=1))

uD<-umap(log(t(D)+1))

pdf(snakemake@output[[1]])
udf<-data.frame(x=uD$layout[,1],y=uD$layout[,2],Cluster=as.factor(clusters[,1]))
p<-ggplot(udf)+geom_point(aes(x,y,color=Cluster))+xlab("UMAP 1")+ylab("UMAP 2")
print(p)
dev.off()
