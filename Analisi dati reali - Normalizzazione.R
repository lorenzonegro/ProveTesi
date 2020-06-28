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

load("Dati Filtrati.RData")
lib.sf.filtered <- librarySizeFactors(filtered)
summary(lib.sf.filtered)

hist(log10(lib.sf.filtered), xlab="Log10[Size factor]", col='grey80')

set.seed(100)
clust.filtered <- quickCluster(filtered) 
table(clust.filtered)

deconv.sf.filtered <- calculateSumFactors(filtered, cluster=clust.filtered)
summary(deconv.sf.filtered)

plot(lib.sf.filtered, deconv.sf.filtered, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(filtered$label)))
abline(a=0, b=1, col="red")

#??? non ho spike ins
#library(scRNAseq)
#filtered <- filtered[,filtered$`single cell quality`=="OK"]
#filtered

#filtered <- computeSpikeFactors(filtered, "ERCC")
#summary(sizeFactors(sce.richard))
#to.plot <- data.frame(
#  DeconvFactor=calculateSumFactors(filtered),
#  SpikeFactor=sizeFactors(filtered),
#  Stimulus=filtered$stimulus, 
#  Time=filtered$time
#)

#ggplot(to.plot, aes(x=DeconvFactor, y=SpikeFactor, color=Time)) +
#  geom_point() + facet_wrap(~Stimulus) + scale_x_log10() + 
#  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")


#Applico fattori di normalizzazione
set.seed(100)
clust.filtered <- quickCluster(filtered) 
filtered <- computeSumFactors(filtered, cluster=clust.filtered, min.mean=0.1)
filtered <- logNormCounts(filtered)
assayNames(filtered)

save(filtered,file="Dati Filtrati.RData")
