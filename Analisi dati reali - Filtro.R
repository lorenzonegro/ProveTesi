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

dati=BaronPancreasData('human')
dati #20125 8569

location <- rowRanges(dati)
is.mito <- any(seqnames(location)=="MT")

library(scater)
df <- perCellQCMetrics(dati, subsets=list(Mito=is.mito))
df

dati <- addPerCellQC(dati, subsets=list(Mito=is.mito))
colnames(colData(dati))

#Identifying low-quality cells With fixed thresholds

#low-quality cells is to apply thresholds on the QC metrics 
#For example, we might consider cells to be low quality if they have library sizes 
#below 100,000 reads;
#express fewer than 5,000 genes; 
#have spike-in proportions above 10%; or have mitochondrial proportions above 10%.

quantile(df$sum,0.05)
qc.lib <- df$sum < 2000 #scarto meno di 2000 reads

quantile(df$detected,0.05)
qc.nexprs <- df$detected < 1000 # scarto meno di 1000 geni

#qc.spike <- df$altexps_ERCC_percent > 10
#qc.mito <- df$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs), Total=sum(discard))
#Scarto 370 cellule per numero di reads, 594 per numero di geni, totale 701


#Identifying outliers
qc.lib2 <- isOutlier(df$sum, log=TRUE, type="lower")
qc.nexprs2 <- isOutlier(df$detected, log=TRUE, type="lower")
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")

#qc.spike2 <- isOutlier(df$altexps_ERCC_percent, type="higher")
#attr(qc.spike2, "thresholds")
#qc.mito2 <- isOutlier(df$subsets_Mito_percent, type="higher")
#attr(qc.mito2, "thresholds")
discard2 <- qc.lib2 | qc.nexprs2

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2),Total=sum(discard2))
#Non scarto per via outliers

#Effetti di batch????
batch <- paste0(dati$donor, "-", dati$Plate)
batch.reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent",
                                                      "altexps_ERCC_percent"), batch=batch)
colSums(as.matrix(batch.reasons))


#metodo alternativo
stats <- cbind(log10(df$sum), log10(df$detected))

library(robustbase)
outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)

#questa riga la scrivo io
reasons.discard <- discard | multi.outlier

#6.5
#Dataset pbmc, come faccio ad usare il mio??
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)
sce.pbmc

#6.6 Removing low-quality cells

# Keeping the columns we DON'T want to discard.
filtered <- dati[,!reasons.discard]

# Using the 'discard' vector for demonstration, 
# as it has more cells for stable calculation of 'lost'.
lost <- calculateAverage(counts(dati)[,!reasons.discard])
kept <- calculateAverage(counts(dati)[,reasons.discard])

library(edgeR)
logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16)

save(filtered,file="Dati Filtrati.RData")
