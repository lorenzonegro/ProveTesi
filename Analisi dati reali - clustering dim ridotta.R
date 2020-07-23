#Analisi dati reali - Clustering
rm(list=ls())
library(devtools)
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/BiocNeighbors")
load_all()
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/clusterExperiment")
load_all()
library(scran)
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/ProveTesi")
library(mclust)
library(igraph)
library(scRNAseq)
library(scater)

load("Dataset finale.RData")
sce.dati <- mnn.out

trueCluster=BaronPancreasData('human')$label
data=t(reducedDim(mnn.out))
cl=clusterMany(data,k=10:20,clusterFunction="kmeans")@clusterMatrix

prop=seq(0,0.9,by=0.1)
adj_CM=apply(cl,2,function(x) adjustedRandIndex(trueCluster,x))

nclust1=st1=idx1=rep(NA,length(prop))
nclust2=st2=idx2=rep(NA,length(prop))
nclust3=st3=idx3=rep(NA,length(prop))
nclust4=st4=idx4=rep(NA,length(prop))
nclust5=st5=idx5=rep(NA,length(prop))

for(i in 1:length(prop))
{
  cat(i)
  #MC1: makeConsensus classico con proportion da 0 a 1
  st1[i] <- system.time(mc1 <- makeConsensus(cl,proportion=prop[i]))[[3]]
  nclust1[i]=length(table(mc1$clustering))
  idx1[i]=adjustedRandIndex(trueCluster,mc1$clustering)
}
nclust1
st1
idx1

for(i in 1:length(prop))
{
  cat(i)
  #MC2: makeConsensus nuovo prop diversi algoritmo walktrap
  st2[i] <- system.time(mc2 <- makeConsensus3(cl, proportion=prop[i], algorithm="cluster_walktrap"))[[3]]
  nclust2[i]=length(table(mc2$clustering))
  idx2[i]=adjustedRandIndex(trueCluster,mc2$clustering)
}
nclust2
idx2
st2

for(i in 1:length(prop))
{
  cat(i)
  #MC3: makeConsensus nuovo prop diversi algoritmo louvain
  st3[i] <- system.time(mc3 <- makeConsensus3(cl, proportion=prop[i], algorithm="cluster_louvain"))[[3]]
  nclust3[i]=length(table(mc3$clustering))
  idx3[i]=adjustedRandIndex(trueCluster,mc3$clustering)
}
nclust3
st3
idx3

for(i in 1:length(prop))
{
  cat(i)
  st4[i] <- system.time(mc4 <- makeConsensus3(cl, proportion=prop[i], algorithm="components_weak"))[[3]]
  nclust4[i]=length(table(mc4$clustering))
  idx4[i]=adjustedRandIndex(trueCluster,mc4$clustering)
}
nclust4
st4
idx4

for(i in 1:length(prop))
{
  cat(i)
  st5[i] <- system.time(mc5 <- makeConsensus3(cl, proportion=prop[i], algorithm="components_strong"))[[3]]
  nclust5[i]=length(table(mc5$clustering))
  idx5[i]=adjustedRandIndex(trueCluster,mc5$clustering)
}
nclust5
st5
idx5


par(mfrow=c(1,3))
plot(nclust1,type="l")
plot(st1,type="l")
plot(idx1,type="l")

plot(nclust2,type="l")
plot(st2,type="l")
plot(idx2,type="l")

plot(nclust3,type="l")
plot(st3,type="l")
plot(idx3,type="l")

par(mfrow=c(1,2))
plot(nclust4,type="l")
plot(st4,type="l")
plot(idx4,type="l")

plot(nclust4,type="l")
plot(st4,type="l")
plot(idx4,type="l")

#confronto 
par(mfrow=c(1,3))
plot(nclust1,col=1,type="l",xlab="proportion",xaxt="n",main="makeConsensus")
axis(1, at=1:length(prop), labels=prop)
plot(nclust2,col=2,type="l",xlab="proportion",xaxt="n",main="makeConsensus3")
axis(1, at=1:length(prop), labels=prop)
points(nclust3,col=3,type="l")
legend("topright",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(nclust4,col=4,type="l",xlab="proportion",xaxt="n",main="makeConsensus3 - components")
axis(1, at=1:length(prop), labels=prop)
points(nclust5,col=6,type="l")
legend("topright",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)

par(mfrow=c(1,3))
plot(idx1,col=1,type="l",xlab="proprtion",xaxt="n",main="makeConsensus")
axis(1, at=1:length(prop), labels=prop)
plot(idx2,col=2,type="l",xlab="proportion",xaxt="n",main="makeConsensus3")
axis(1, at=1:length(prop), labels=prop)
points(idx3,col=3,type="l")
legend("bottomright",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(idx4,col=4,type="l",xlab="proportion",xaxt="n",main="makeConsensus3")
axis(1, at=1:length(prop), labels=prop)
points(idx5,col=6,type="l")
legend("bottomright",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)

par(mfrow=c(1,3))
plot(st1,col=1,type="l",xlab="proprtion",xaxt="n",main="makeConsensus")
axis(1, at=1:length(prop), labels=prop)
plot(st2,col=2,type="l",xlab="proportion",xaxt="n",main="makeConsensus3")
axis(1, at=1:length(prop), labels=prop)
points(st3,col=3,type="l")
legend("topright",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(st4,col=4,type="l",xlab="proportion",xaxt="n",main="makeConsensus3")
axis(1, at=1:length(prop), labels=prop)
points(st5,col=6,type="l")
legend("topright",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)

#Confronto migliori prestazioni
#guardo system time quando l'indice è 1
c(max(idx1),prop[which.max(idx1)],st1[which.max(idx1)]) #mc1
c(max(idx2),prop[which.max(idx2)],st2[which.max(idx2)]) #mc3 walktrap
c(max(idx3),prop[which.max(idx3)],st3[which.max(idx3)]) #mc3 louvain
c(max(idx4),prop[which.max(idx4)],st4[which.max(idx4)]) #mc3 c weak
c(max(idx5),prop[which.max(idx5)],st5[which.max(idx5)]) #mc3 c strong

