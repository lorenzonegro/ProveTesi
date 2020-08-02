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
set.seed(11)
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

# 0.74 10-25 0.3 11

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
riassunto=0
for(i in 1:ncol(tab.ridotta))
{
  j=which.max(tabella.class[,i])
  tabella.finale.ridotta[,j]=tabella.finale.ridotta[,j]+tabella.class[,i]
  riassunto=c(riassunto,j)
} 
sum(diag(tabella.finale.ridotta[,1:14]))/sum(tabella.finale)
rownames(tabella.finale.ridotta)=row.names(tabella.class)
write.table(tabella.finale.ridotta,file="Tabella ridotta.txt")


clustering.finale=rep(0,length(trueCluster))
riassunto=riassunto[-1]
riassunto
table(riassunto)
outliers=sce.dati[,colSums(tabella.class)<1]
getBestFeatures(mc3)
tabella.class[,1]
tabella.class[,3]

cluster.aggregati=rep(NA,length(trueCluster))
cluster.aggregati=mc3$clustering

library(plyr)
cluster.aggregati=as.character(cluster.aggregati)
cluster.aggregati <- revalue(cluster.aggregati,c("1"=names(riassunto)[1],"2"=names(riassunto)[2],
                                                 "3"=names(riassunto)[3],"4"=names(riassunto)[4],
                                                 "5"=names(riassunto)[5],"6"=names(riassunto)[6],
                                                 "7"=names(riassunto)[7],"8"=names(riassunto)[8],
                                                 "9"=names(riassunto)[9],"10"=names(riassunto)[10],
                                                 "11"=names(riassunto)[11],"12"=names(riassunto)[12],
                                                 "13"=names(riassunto)[13],"14"=names(riassunto)[14],
                                                 "15"=names(riassunto)[15],"16"=names(riassunto)[16],
                                                 "17"=names(riassunto)[17],"18"=names(riassunto)[18],
                                                 "19"=names(riassunto)[19],"20"=names(riassunto)[20],
                                                 "21"=names(riassunto)[21],"22"=names(riassunto)[22],
                                                 "23"=names(riassunto)[23],"24"=names(riassunto)[24],
                                                 "25"=names(riassunto)[25],"26"=names(riassunto)[26],
                                                 "27"=names(riassunto)[27],"28"=names(riassunto)[28]))
                                         

for(i in 1:length(trueCluster))
{
  if(mc3$clustering[i]>28)
    cluster.aggregati[i]=-1
}
table(cluster.aggregati)
table(trueCluster)

#cluster.aggregati <- revalue(cluster.aggregati,c("acinar"=1,"activated_stellate"=2,"beta"=3,
#                                                 "delta"=4,"ductal"=5,"endothelial"=6,
#                                                 "gamma"=7,"macrophage"=8))
                                                
#per il getbestfeatures devo usare dati originali
load("Dati normalizzati.RData")
provadati=as.matrix(dati@assays@data@listData[["logcounts"]])
mark=findMarkers(provadati,groups=cluster.aggregati)
table(cluster.aggregati)
mark

nomi.gruppi=unique(cluster.aggregati)[-c(7,11)]
nomi.gruppi
geniscelti=NA
for(i in nomi.gruppi)
{
  cat(i)
  chosen <- i
  interesting <- mark[[chosen]]
  geniscelti=c(geniscelti,rownames(interesting[1:10,1:4]))
}
geniscelti=geniscelti[-1]
geniscelti
which(is.na(geniscelti))
geniscelti=unique(na.omit(geniscelti))
tab.heatmap=matrix(NA,length(unique(nomi.gruppi)),length(geniscelti))
dim(tab.heatmap)
j=0
for(i in nomi.gruppi)
{
 j=j+1
 tab.heatmap[j,]=rowMeans(provadati[geniscelti,which(cluster.aggregati==i)])
} 
rownames(tab.heatmap)=nomi.gruppi
colnames(tab.heatmap)=geniscelti

?pheatmap
pheatmap(tab.heatmap)
