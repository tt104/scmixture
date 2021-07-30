logsum<-function(a,b)
{
	m<-pmax(a,b)
	log(exp(a-m)+exp(b-m))+m
}

maxlikeNB2<-function(xs,scales)
{
	mu<-NULL
	os<-NULL
	ws<-NULL
	
	for(j in 1:ncol(xs))
	{
		mask<-which(xs[,j]>0)
		zmask<-which(xs[,j]==0)
		nblikeZI<-function(ps)
		{
			if(length(zmask>0))
			{
				l0<- -sum(logsum(log(ps[3])+dnbinom(xs[zmask,j],size=ps[1],mu=scales[zmask]*ps[2],log=TRUE),rep(log(1-ps[3]),length(zmask))))
			}
			else
			{
				l0<-0.0
			}
			l<- -sum(log(ps[3])+dnbinom(xs[mask,j],size=ps[1],mu=scales[mask]*ps[2],log=TRUE))
			p<- -dgamma(ps[1],shape=1,rate=0.01,log=TRUE)
			l0+l+p
		}
		mlzi<-optim(c(1,mean(xs[,j]),0.5),nblikeZI,method="L-BFGS-B",lower=c(0.001,0.001,0.001),upper=c(Inf,Inf,0.9999))
		mu<-c(mu,mlzi$par[2])
		os<-c(os,mlzi$par[1])
		ws<-c(ws,mlzi$par[3])
	}
	list(mu=mu,os=os,ws=ws)
}

d<-as.matrix(read.csv(snakemake@input[[1]],header=TRUE,row.names=1))

dd<-t(d)

ldata<-log(dd+1)
vs<-apply(ldata,2,var)
cv<-sqrt(exp(vs-1))
ix<-order(cv,decreasing=TRUE)[1:1000]

xs<-t(d)[,ix[1:1000]]

scales<-as.matrix(read.csv(snakemake@input[[2]],header=FALSE))

res<-maxlikeNB2(xs,scales)

ms<-res$mu
os<-res$os
ws<-res$ws

gmulike<-function(ps)
{
	        -sum(dgamma(ms,shape=ps[1],rate=ps[2],log=TRUE))
}
             
mug<-optim(c(mean(ms),1),gmulike,method="L-BFGS-B",lower=c(0.00001,0.00001),upper=c(Inf,Inf))

gomegalike<-function(ps)
{
	        -sum(dgamma(os,shape=ps[1],rate=ps[2],log=TRUE))
}
             
omegag<-optim(c(mean(os),1),gomegalike,method="L-BFGS-B",lower=c(0.00001,0.00001),upper=c(Inf,Inf))

hp<-c(mug$par,omegag$par)

write.table(hp,snakemake@output[[1]],row.names=FALSE,col.names=FALSE)
