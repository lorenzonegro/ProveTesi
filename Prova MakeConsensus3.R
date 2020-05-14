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

#Per usare adjuster rand index devo dire ogni oss a che cluster appartiene

load("simCount_10.RDAta")
load("simData_10.RData")
load("trueCluster_10.RData")


cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix

mc3=makeConsensus3(cl,proportion=0.1,"cluster_walktrap") #ok ci sta ma dovrebbe darmi errore

table(mc3$clustering)

mc3=makeConsensus3(cl,proportion=0.5,"cluster_walktrap") #ok ci sta ma dovrebbe darmi errore

table(mc3$clustering)

mc3=makeConsensus3(cl,proportion=0.7,"cluster_walktrap") #ok ci sta ma dovrebbe darmi errore

table(mc3$clustering)

mc3=makeConsensus3(cl,proportion=1,"cluster_walktrap") #ok ci sta ma dovrebbe darmi errore

table(mc3$clustering)
