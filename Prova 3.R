rm(list=ls())
library(devtools)

setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/BiocNeighbors")
load_all()

setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/clusterExperiment")
load_all()

data=matrix(c(1,1,1,
              1,1,2,
              3,1,2,
              2,2,2,
              1,1,1),5,3,byrow=T)

library(scran)
buildSNNGraph(data, BNPARAM = VptreeParam(distance = "Hamming"),transposed=TRUE)

cl=read.csv("example_cl.csv")
head(cl)
g=buildSNNGraph(cl,k=5, BNPARAM = VptreeParam(distance = "Hamming"),transposed=TRUE)
#RICORDA: scelta di d
dim(cl)
plot(g)

library(igraph)
cluster_g=cluster_walktrap(g)
names(cluster_g)
cluster_g$membership

#library(clusterExperiment)
cl=as.matrix(cl)
mc1=makeConsensus(cl,proportion=0.6)

mc2=makeConsensus2(cl,k=50,algorithm="cluster_louvain")

table(mc1$clustering)
table(mc2$clustering)

which.max(mc1$clustering)
which.max(mc2$clustering)
#aggiungi system.time

system.time(makeConsensus(cl,proportion=1))
            
system.time(makeConsensus2(cl,k=10,algorithm="cluster_walktrap"))


