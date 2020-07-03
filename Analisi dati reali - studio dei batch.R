#Analisi dati reali - Riduzione dimensionalità
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

load("Dati con PCA.RData")

# Using RandomParam() as it is more efficient for file-backed matrices.
set.seed(0010101010)
uncorrected <- runPCA(sce.dati,BSPARAM=BiocSingular::RandomParam())
snn.gr <- buildSNNGraph(uncorrected, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters, Batch=uncorrected$donor)
tab
set.seed(1111001)
uncorrected <- runTSNE(uncorrected, dimred="PCA")
plotTSNE(uncorrected, colour_by="donor")

set.seed(1000101001)
library(batchelor)
mnn.out <- fastMNN(sce.dati$donor, d=50, k=20,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out