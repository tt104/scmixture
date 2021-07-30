# Original data
D<-as.matrix(read.csv(snakemake@input[[1]],header=TRUE,row.names=1))
cellnames<-colnames(D)

# Gene names selected
names<-as.matrix(read.csv(snakemake@input[[5]],header=FALSE))[,1]
ngene<-length(names)

# Re-map clusters
clustersRaw<-as.matrix(read.csv(snakemake@input[[2]],header=FALSE))
clist<-unique(clustersRaw)

nc<-length(clist)

clusters<-rep(0,length(clustersRaw))
clustersDict<-list()

for(i in 1:nc)
{
	t<-clist[i]
	clusters[clustersRaw==t]<-i
	clustersDict[[t]]<-i
}

# Write counts table
write.table(clusters,snakemake@output[[1]],sep=',',row.names=cellnames,col.names=c("Cluster"),quote=FALSE)
write.table(table(clusters),snakemake@output[[2]],sep=',')

# Fetch cluster means and dispersions
clusterMuRaw<-as.matrix(read.csv(snakemake@input[[3]],header=FALSE))
clusterOmegaRaw<-as.matrix(read.csv(snakemake@input[[4]],header=FALSE))


clusterMu<-matrix(0,nrow=nc,ncol=ngene)
colnames(clusterMu)<-names
rownames(clusterMu)<-c(1:nc)

clusterOmega<-matrix(0,nrow=nc,ncol=ngene)
colnames(clusterOmega)<-names
rownames(clusterOmega)<-c(1:nc)

for(i in 1:nc)
{
	t<-clist[i]
	s<-clustersDict[[t]]
	clusterMu[s,]<-clusterMuRaw[t,]
	clusterOmega[s,]<-clusterOmegaRaw[t,]
}

write.table(clusterMu,snakemake@output[[3]],sep=',',quote=FALSE)
write.table(clusterOmega,snakemake@output[[4]],sep=',',quote=FALSE)
