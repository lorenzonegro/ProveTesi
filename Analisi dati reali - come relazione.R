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
library(GEOquery)
library(recount)
library(org.Hs.eg.db)
library(EDASeq)
library(gridExtra)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(clusterProfiler)
library(globaltest)

dati=BaronPancreasData('human')
dati #20125 8569

# Creiamo il summarized experiment
se <- SummarizedExperiment(assays = list(counts=dati@assays@data@listData[["counts"]]),
colData = dati@colData)
rse_gene=se

ids <- mapIds(org.Hs.eg.db, rownames(se), "ENTREZID", "SYMBOL")
rownames(se)<-ids
se1<-se[!is.na(ids),]

# filtraggio dei geni poco espressi
filter <-rowMeans(assay(se1))>=0.001
filtered <- se1[filter,]

pal <-brewer.pal(8,"Set1")

assay(filtered,"uq")<-betweenLaneNormalization(assay(filtered),which="upper")
assay(filtered,"fq")<-betweenLaneNormalization(assay(filtered),which="full")

