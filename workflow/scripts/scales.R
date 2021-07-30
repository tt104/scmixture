library(scran)

D<-as.matrix(read.csv(snakemake@input[[1]],header=TRUE,row.names=1))

clusters <- quickCluster(D)
sce <- calculateSumFactors(D, clusters=clusters)

write.table(sce,snakemake@output[[1]],sep=',',row.names=FALSE,col.names=FALSE)
