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

analisi <- function(ngruppi=10,range=c(-155.5,155.5),n_cells=3000,n_genes=150)
{
  medie_x=runif(ngruppi,range[1],range[2])
  medie_y=rev(medie_x)
  prop_gruppi=c(runif(ngruppi,0,1))
  prop_gruppi=prop_gruppi/sum(prop_gruppi)
  
  simulazioni=simulate_gauss_mix(n_cells=n_cells, n_genes=n_genes,
                                 k=ngruppi, x_mus = medie_x, 
                                 x_sds = c(runif(ngruppi,0.1,1)), 
                                 y_mus = medie_y, 
                                 y_sds = c(runif(ngruppi,0.1,1)), 
                                 prop1 = prop_gruppi)
  
  trueCluster=simulazioni$true_cluster_id
  simData=as.matrix(t(simulazioni$obs_data))
  
  prop=0.5
  
  nclust1=st1=idx1=rep(NA,length(prop))
  nclust3=st3=idx3=rep(NA,length(prop))
  
  cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix
  
  #MC1: makeConsensus classico con proportion da 0 a 1
  st1 <- system.time(mc1 <- makeConsensus(cl,proportion=0.5))[[3]]
  nclust1=length(table(mc1$clustering))
  idx1=adjustedRandIndex(trueCluster,mc1$clustering)
  
  #MC3: makeConsensus nuovo prop diversi algoritmo louvain
  st3 <- system.time(mc3 <- makeConsensus3(cl, proportion=0.5, algorithm="cluster_louvain"))[[3]]
  nclust3=length(table(mc3$clustering))
  idx3=adjustedRandIndex(trueCluster,mc3$clustering)
  
  #restituisco in output una matrice con numero di righe pari a proportion (10)
  #3 colonne, una per system time, una per nclust, una per l'indice
  #questa matrice è una per ognuno dei metodi, quindi 5 in totale
  mc1=cbind(nclust1,st1,idx1)
  mc_louvain=cbind(nclust3,st3,idx3)
  return(list(mc=mc1,mc_louvain=mc_louvain))
}

sim=analisi(ngruppi=10,n_genes = 50,range=c(-155.5,155.5))
nclust_mc=sim[[1]][,1]
nclust_lo=sim[[2]][,1]

st_mc=sim[[1]][,2]
st_lo=sim[[2]][,2]

idx_mc=sim[[1]][,3]
idx_lo=sim[[2]][,3]

numero_geni=c(150,300,500,1000)

#prima simulazione casuale
set.seed(456)
for(i in 1:length(numero_geni))
{
  cat(i)
  sim=analisi(ngruppi=10,n_cells=3000,n_genes=numero_geni[i],range=c(-155.5,155.5))
  nclust_mc=cbind(nclust_mc,sim[[1]][,1])
  nclust_lo=cbind(nclust_lo,sim[[2]][,1])
  
  st_mc=cbind(st_mc,sim[[1]][,2])
  st_lo=cbind(st_lo,sim[[2]][,2])
  
  idx_mc=cbind(idx_mc,sim[[1]][,3])
  idx_lo=cbind(idx_lo,sim[[2]][,3])
}

output=list("Numero di cluster mC"=nclust_mc,
            "Numero di cluster Louvain"=nclust_lo,
            "System time mC"=st_mc,
            "System time Louvain"=st_lo,
            "Indice mc"=idx_mc,
            "Indice Louvain"=idx_lo)

save(output,file="Statistiche cambio numero geni.RData")


