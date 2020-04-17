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

cl=clusterMany(simData,k=2:10,clusterFunction="kmeans")@clusterMatrix #k numero cluster
head(cl)
g=buildSNNGraph(cl,k=10, BNPARAM = VptreeParam(distance = "Hamming"),transposed=T) #k di buildgraph numero di vicini + vicini
#RICORDA: scelta di d
dim(cl)
plot(g)

library(igraph)
cluster_g=cluster_walktrap(g)
names(cluster_g)
cluster_g$membership

#library(clusterExperiment)
mc1=makeConsensus(cl,proportion=0.7)

mc2=makeConsensus2(cl,k=10,algorithm="cluster_walktrap") #aggiungere step a makeconsensus2

table(mc1$clustering)
table(mc2$clustering)
trueCluster
#aggiungi system.time

system.time(makeConsensus(cl,proportion=0.2))

system.time(makeConsensus2(cl,k=15,algorithm="cluster_louvain"))

#Quanto cambia se cambio make consensus? Cluster many sembra molto più importante
#make consensus sceglie il migliore?
#Sembra scegliere quello con più cluster

#Mergecluster x prossima volta
#Valutare anche system time e uso della RAM
#Rprof x utilizzo ram
#simulare da normale trivariata in stile ?simData 3000 dati al posto di 300
#adjusted rand index nel pachetto mclust anche se il numero di cluster non corrisponde