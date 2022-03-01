library(mclust)

labels<-factor(as.vector(as.matrix(read.csv(args[1],header=FALSE))))

z<-as.matrix(read.csv(snakemake@input[[2]],header=TRUE))[,1]

ari<-adjustedRandIndex(as.integer(z),as.integer(labels))
res<-data.frame(ari=ari)
write.table(res,snakemake@output[[1]])
