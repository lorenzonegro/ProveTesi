rm(list=ls())

load("Statistiche 3 gruppi.RData")

nclust_mc=round(rowMeans(output$`Numero di cluster mC`))
nclust_wa=round(rowMeans(output$`Numero di cluster walktrap`))
nclust_lo=round(rowMeans(output$`Numero di cluster Louvain`))
nclust_cw=round(rowMeans(output$`Numero di cluster walktrap`))
nclust_cs=round(rowMeans(output$`Numero di cluster comp_strong`))

prop=seq(0,0.9,by=0.1)

#confronto 
par(mfrow=c(1,3))
plot(nclust_mc,col=1,type="l",xlab="proportion",xaxt="n",main="makeConsensus")
axis(1, at=1:length(prop), labels=prop)
plot(nclust_wa,col=2,type="l",xlab="proportion",xaxt="n",main="makeConsensus3 - Louvain & walktrap")
axis(1, at=1:length(prop), labels=prop)
points(nclust_lo,col=3,type="l")
legend("topleft",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(nclust_cw,col=4,type="l",xlab="proportion",xaxt="n",main="makeConsensus3 - components")
axis(1, at=1:length(prop), labels=prop)
points(nclust_cs,col=6,type="l")
legend("topleft",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)




par(mfrow=c(1,2))
plot(idx1,col=1,type="l",xlab="proprtion",xaxt="n",main="makeConsensus")
axis(1, at=1:length(prop), labels=prop)
plot(idx2,col=2,type="l",xlab="proportion",xaxt="n",main="makeConsensus2")
axis(1, at=1:length(prop), labels=prop)
points(idx3,col=3,type="l")
legend("topleft",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(idx4,col=4,type="l",xlab="proportion",xaxt="n",main="makeConsensus2")
axis(1, at=1:length(prop), labels=prop)
points(idx5,col=6,type="l")
legend("bottomleft",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)
plot(adj_CM,col=4,type="l",xlab="k",xaxt="n",main="risultato cluster many")
axis(1, at=1:19, labels=c(2:20))


par(mfrow=c(1,3))
plot(st1,col=1,type="l",xlab="proprtion",xaxt="n",main="makeConsensus")
axis(1, at=1:length(prop), labels=prop)
plot(st2,col=2,type="l",xlab="proportion",xaxt="n",main="makeConsensus2")
axis(1, at=1:length(prop), labels=prop)
points(st3,col=3,type="l")
legend("topleft",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(st4,col=4,type="l",xlab="proportion",xaxt="n",main="makeConsensus2")
axis(1, at=1:length(prop), labels=prop)
points(st5,col=6,type="l")
legend("topleft",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)