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
cl=clusterMany(data,k=10:25,clusterFunction="kmeans")@clusterMatrix

prop=0.3
adj_CM=apply(cl,2,function(x) adjustedRandIndex(trueCluster,x))

nclust3=st3=idx3=rep(NA,length(prop))

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

# 0.74 10-25 0.3

tabella.class=as.matrix(table(trueCluster,mc3$clustering))
dim(tabella.class)
write.table(tabella.class,file="Tabella iniziale.txt")

length(table(trueCluster))
tabella.finale=matrix(0,length(table(trueCluster)),length(table(trueCluster)))
colnames(tabella.finale)=row.names(tabella.class)
rownames(tabella.finale)=row.names(tabella.class)
for(i in 1:ncol(tabella.class))
{
  j=which.max(tabella.class[,i])
  tabella.finale[,j]=tabella.finale[,j]+tabella.class[,i]
}  

sum(diag(tabella.finale))/sum(tabella.finale)

#alternativa
tab.ridotta=tabella.class[,(colSums(tabella.class)>1)]
tabella.finale.ridotta=matrix(0,dim(tab.ridotta)[1],dim(tab.ridotta)[2])
for(i in 1:ncol(tab.ridotta))
{
  j=which.max(tabella.class[,i])
  tabella.finale.ridotta[,j]=tabella.finale.ridotta[,j]+tabella.class[,i]
} 
sum(diag(tabella.finale.ridotta[,1:14]))/sum(tabella.finale)
rownames(tabella.finale.ridotta)=row.names(tabella.class)
write.table(tabella.finale.ridotta,file="Tabella ridotta.txt")
