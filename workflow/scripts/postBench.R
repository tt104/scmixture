# Re-map clusters
clustersRaw<-as.matrix(read.csv(snakemake@input[[1]],header=FALSE))
clist<-unique(clustersRaw)

nc<-length(clist)

clusters<-rep(0,length(clustersRaw))

for(i in 1:nc)
{
	t<-clist[i]
	clusters[clustersRaw==t]<-i
}

# Write counts table
write.table(clusters,snakemake@output[[1]],sep=',',row.names=cellnames,col.names=c("Cluster"),quote=FALSE)
