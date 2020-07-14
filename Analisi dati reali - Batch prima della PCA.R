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

load("Dati con geni più significativi.RData")
sce.dati <- sce.dati.hvg
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
mnn.out <- fastMNN(sce.dati, batch=sce.dati$donor, d=50, k=20, BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out
dim(reducedDim(mnn.out, "corrected"))

snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
tab.mnn

set.seed(0010101010)
mnn.out <- runTSNE(mnn.out, dimred="corrected")

mnn.out$batch <- factor(mnn.out$batch)
plotTSNE(mnn.out, colour_by="batch")

m.out <- findMarkers(uncorrected, clusters.mnn, block=uncorrected$batch,
                     direction="up", lfc=1, row.data=rowData(uncorrected)[,3,drop=FALSE])

save(mnn.out,file="Dataset batch prima della PCA.RData")
