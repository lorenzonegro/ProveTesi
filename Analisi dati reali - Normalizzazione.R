#Analisi dati reali - Normalizzazioni
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

dati=BaronPancreasData('human')
lib.sf.dati <- librarySizeFactors(dati)
summary(lib.sf.dati)

hist(log10(lib.sf.dati), xlab="Log10[Size factor]", col='grey80')

set.seed(100)
clust.dati <- quickCluster(dati) 
table(clust.dati)

deconv.sf.dati <- calculateSumFactors(dati, cluster=clust.dati)
summary(deconv.sf.dati)

plot(lib.sf.dati, deconv.sf.dati, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(dati$label)))
abline(a=0, b=1, col="red")

#non ho spike ins

#Applico fattori di normalizzazione
set.seed(100)
clust.dati <- quickCluster(dati) 
dati <- computeSumFactors(dati, cluster=clust.dati, min.mean=0.1)
dati <- logNormCounts(dati)
assayNames(dati)

save(dati,file="Dati Normalizzati.RData")

library(scran)
dec.dati <- modelGeneVar(dati)

# Visualizing the fit:
fit.dati <- metadata(dec.dati)
plot(fit.dati$mean, fit.dati$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.dati$trend(x), col="dodgerblue", add=TRUE, lwd=2)

dec.dati[order(dec.dati$bio, decreasing=TRUE),] 

hvg.dati.var <- getTopHVGs(dec.dati, n=1000)
str(hvg.dati.var)
hvg.dati.var.2 <- getTopHVGs(dec.dati, fdr.threshold=0.05)
length(hvg.dati.var.2)

#scelgo secondo metodo
chosen <- hvg.dati.var.2
sce.dati.hvg <- dati[chosen,]
dim(sce.dati.hvg)

save(sce.dati.hvg,file="Dati con geni più significativi")
