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

simulate_gauss_mix <- function(n_cells, n_genes,
                               k, x_mus = c(0,5,5), 
                               x_sds = c(1,0.1,1), 
                               y_mus = c(5,5,0), 
                               y_sds = c(1,0.1,1), 
                               prop1 = c(0.3,0.5,0.2))
{ 
  
  if(k != length(x_mus)){stop("k is not same as length of x_mus")} 
  if(k != length(x_sds)){stop("k is not same as length of x_sds")} 
  if(k != length(y_mus)){stop("k is not same as length of y_mus")} 
  if(k != length(y_sds)){stop("k is not same as length of y_sds")} 
  if(k != length(prop1)){stop("k is not same as length of prop1")} 
  
  comp1 <- sample(seq_len(k), prob=prop1, size=n_cells, replace=TRUE)
  
  # Sampling locations for cells in each component
  samples1 <- cbind(rnorm(n=n_cells, mean=x_mus[comp1],sd=x_sds[comp1]),
                    rnorm(n=n_cells, mean=y_mus[comp1],sd=y_sds[comp1]))
  
  # Random projection to D dimensional space, to mimic high-dimensional expression data.
  proj <- matrix(rnorm(n_genes*n_cells), nrow=n_genes, ncol=2)
  A1 <- samples1 %*% t(proj)
  
  # Add normally distributed noise.
  A1 <- A1 + rnorm(n_genes*n_cells)
  rownames(A1) <- paste0("Cell", seq_len(n_cells), "-1")
  colnames(A1) <- paste0("Gene", seq_len(n_genes))
  
  list("true_center" = cbind("x" = x_mus, "y" = y_mus),
       "true_cluster_id" = comp1,
       "true_data" = samples1, 
       "obs_data" = A1)
}

ngruppi=10
medie_x=runif(ngruppi,-155.5,155.5)
medie_y=rev(medie_x)
prop_gruppi=c(runif(10,0,1))
prop_gruppi=prop_gruppi/sum(prop_gruppi)

simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=50,
                               k=ngruppi, x_mus = medie_x, 
                               x_sds = c(runif(ngruppi,0.1,1)), 
                               y_mus = medie_y, 
                               y_sds = c(runif(ngruppi,0.1,1)), 
                               prop1 = prop_gruppi)

trueCluster=simulazioni$true_cluster_id
simData=as.matrix(t(simulazioni$obs_data))

cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix

Rprof(filename ="ram_consensus.txt", append = FALSE, memory.profiling = TRUE)
mc1=makeConsensus(cl,proportion=0.5)
Rprof(NULL)

profile <- summaryRprof(filename = "ram_consensus.txt", chunksize = -1L, 
                        memory = "tseries", diff = FALSE)
max_mem_cons1 <- max(rowSums(profile[,1:3]))*0.00000095367432
consumo_memoria=max_mem_cons1

#150
simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=150,
                               k=ngruppi, x_mus = medie_x, 
                               x_sds = c(runif(ngruppi,0.1,1)), 
                               y_mus = medie_y, 
                               y_sds = c(runif(ngruppi,0.1,1)), 
                               prop1 = prop_gruppi)

trueCluster=simulazioni$true_cluster_id
simData=as.matrix(t(simulazioni$obs_data))

cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix

Rprof(filename ="ram_consensus.txt", append = FALSE, memory.profiling = TRUE)
mc1=makeConsensus(cl,proportion=0.5)
Rprof(NULL)

profile <- summaryRprof(filename = "ram_consensus.txt", chunksize = -1L, 
                        memory = "tseries", diff = FALSE)
max_mem_cons1 <- max(rowSums(profile[,1:3]))*0.00000095367432
consumo_memoria=c(consumo_memoria,max_mem_cons1)

#300
simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=300,
                               k=ngruppi, x_mus = medie_x, 
                               x_sds = c(runif(ngruppi,0.1,1)), 
                               y_mus = medie_y, 
                               y_sds = c(runif(ngruppi,0.1,1)), 
                               prop1 = prop_gruppi)

trueCluster=simulazioni$true_cluster_id
simData=as.matrix(t(simulazioni$obs_data))

cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix

Rprof(filename ="ram_consensus.txt", append = FALSE, memory.profiling = TRUE)
mc1=makeConsensus(cl,proportion=0.5)
Rprof(NULL)

profile <- summaryRprof(filename = "ram_consensus.txt", chunksize = -1L, 
                        memory = "tseries", diff = FALSE)
max_mem_cons1 <- max(rowSums(profile[,1:3]))*0.00000095367432
consumo_memoria=c(consumo_memoria,max_mem_cons1)

#500
simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=500,
                               k=ngruppi, x_mus = medie_x, 
                               x_sds = c(runif(ngruppi,0.1,1)), 
                               y_mus = medie_y, 
                               y_sds = c(runif(ngruppi,0.1,1)), 
                               prop1 = prop_gruppi)

trueCluster=simulazioni$true_cluster_id
simData=as.matrix(t(simulazioni$obs_data))

cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix

Rprof(filename ="ram_consensus.txt", append = FALSE, memory.profiling = TRUE)
mc1=makeConsensus(cl,proportion=0.5)
Rprof(NULL)

profile <- summaryRprof(filename = "ram_consensus.txt", chunksize = -1L, 
                        memory = "tseries", diff = FALSE)
max_mem_cons1 <- max(rowSums(profile[,1:3]))*0.00000095367432
consumo_memoria=c(consumo_memoria,max_mem_cons1)

#1000
simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=1000,
                               k=ngruppi, x_mus = medie_x, 
                               x_sds = c(runif(ngruppi,0.1,1)), 
                               y_mus = medie_y, 
                               y_sds = c(runif(ngruppi,0.1,1)), 
                               prop1 = prop_gruppi)

trueCluster=simulazioni$true_cluster_id
simData=as.matrix(t(simulazioni$obs_data))

cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix

Rprof(filename ="ram_consensus.txt", append = FALSE, memory.profiling = TRUE)
mc1=makeConsensus(cl,proportion=0.5)
Rprof(NULL)

profile <- summaryRprof(filename = "ram_consensus.txt", chunksize = -1L, 
                        memory = "tseries", diff = FALSE)
max_mem_cons1 <- max(rowSums(profile[,1:3]))*0.00000095367432
consumo_memoria=c(consumo_memoria,max_mem_cons1)

load("Statistiche cambio numero geni.RData")
output[["Memoria mC"]] <- consumo_memoria
save(output,file="Statistiche cambio numero geni.RData")