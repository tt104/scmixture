dir.create(snakemake@output[[1]])

# Fetch cluster means and dispersions
clusterMu<-as.matrix(read.csv(snakemake@input[[1]],header=TRUE,row.names=1))
clusterOmega<-as.matrix(read.csv(snakemake@input[[2]],header=TRUE,row.names=1))

nc<-nrow(clusterMu)
names<-colnames(clusterMu)

avmu<-colMeans(clusterMu)
avom<-colMeans(clusterOmega)
qup<-qnbinom(0.95,size=avom,mu=avmu)
qlo<-qnbinom(0.05,size=avom,mu=avmu)

for(cl in c(1:nc))
{
	clmu<-clusterMu[cl,]
	clom<-clusterOmega[cl,]

	up<-which(clmu>qup)
	down<-which(clmu<qlo)
	upgenes<-names[up]
	downgenes<-names[down]

	if(length(upgenes)!=0)
	{
	write.table(upgenes,paste(snakemake@output[[1]],"/Cluster_",cl,"_up.csv",sep=''),sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
	if(length(downgenes)!=0)
	{
	write.table(downgenes,paste(snakemake@output[[1]],"/Cluster_",cl,"_down.csv",sep=''),sep=',',row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
}
