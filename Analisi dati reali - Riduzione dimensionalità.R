#Analisi dati reali - Riduzione dimensionalit�
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

load("Dati con geni pi� significativi.RData")
dim(sce.dati.hvg)

set.seed(100) # See below.
sce.dati <- runPCA(sce.dati.hvg) 
reducedDimNames(sce.dati)

dim(reducedDim(sce.dati, "PCA"))

# Percentage of variance explained is tucked away in the attributes.
percent.var <- attr(reducedDim(sce.dati), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")

reducedDim(sce.dati, "PCA") <- reducedDim(sce.dati, "PCA")[,1:48]
ncol(reducedDim(sce.dati, "PCA"))

dim(sce.dati)
dim(reducedDim(sce.dati,"PCA"))

