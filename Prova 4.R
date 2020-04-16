#Prova4 
#Test con SimData utilizzando anche clusterMany
rm(list=ls())
library(devtools)

setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/BiocNeighbors")
load_all()

setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/clusterExperiment")
load_all()

library(scran)
data("simData")

cl=clusterMany(simData,k=2:3,clusterFunction="kmeans")@clusterMatrix
head(cl)
g=buildSNNGraph(cl,k=10, BNPARAM = VptreeParam(distance = "Hamming"),transposed=T)
#RICORDA: scelta di d
dim(cl)
plot(g)

library(igraph)
cluster_g=cluster_walktrap(g)
names(cluster_g)
cluster_g$membership

#library(clusterExperiment)
mc1=makeConsensus(cl,proportion=0.9)

mc2=makeConsensus2(cl,k=10,algorithm="cluster_louvain")

table(mc1$clustering)
table(mc2$clustering)
trueCluster
#aggiungi system.time

system.time(makeConsensus(cl,proportion=1))

system.time(makeConsensus2(cl,k=10,algorithm="cluster_walktrap"))
