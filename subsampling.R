rm(list=ls())
library(devtools)
#install_github("lorenzonegro/BiocNeighbors")
#install_github("epurdom/clusterExperiment@biocneighbors")
library(BiocNeighbors)
library(clusterExperiment)
library(scran)
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
range=c(-10.5,10.5)
medie_x=runif(ngruppi,range[1],range[2])
medie_y=rev(medie_x)
prop_gruppi=c(runif(ngruppi,0,1))
prop_gruppi=prop_gruppi/sum(prop_gruppi)
n_cells=3000
n_genes=150

simulazioni=simulate_gauss_mix(n_cells=n_cells, n_genes=n_genes,
                               k=ngruppi, x_mus = medie_x, 
                               x_sds = c(runif(ngruppi,0.1,1)), 
                               y_mus = medie_y, 
                               y_sds = c(runif(ngruppi,0.1,1)), 
                               prop1 = prop_gruppi)

trueCluster=simulazioni$true_cluster_id
simData=as.matrix(t(simulazioni$obs_data))

ce=clusterSingle(simData, subsample = TRUE, sequential = FALSE,
              mainClusterArgs = list(clusterFunction = "hierarchical01",
                                     clusterArgs = list(alpha = 1)),
              subsampleArgs = list(clusterFunction = "kmeans",
                                   clusterArgs = list(k = 5),
                                   samp.p = 0.7,
                                   resamp.num = 100))

table(ce@clusterMatrix)

nclust1=st1=idx1=matrix(NA,10,100)
alpha.val=seq(0.1:0.9,by=0.1)
analisi <- function(ngruppi,range,n_cells=3000,n_genes=150)
{
  medie_x=runif(ngruppi,range[1],range[2])
  medie_y=rev(medie_x)
  prop_gruppi=c(runif(ngruppi,0,1))
  prop_gruppi=prop_gruppi/sum(prop_gruppi)
  n_cells=3000
  n_genes=150
  simulazioni=simulate_gauss_mix(n_cells=n_cells, n_genes=n_genes,
                                 k=ngruppi, x_mus = medie_x, 
                                 x_sds = c(runif(ngruppi,0.1,1)), 
                                 y_mus = medie_y, 
                                 y_sds = c(runif(ngruppi,0.1,1)), 
                                 prop1 = prop_gruppi)
  
  trueCluster=simulazioni$true_cluster_id
  simData=as.matrix(t(simulazioni$obs_data))
  
  prop=seq(0,0.9,by=0.1)
  nclust1=st1=idx1=rep(NA,length(prop))
  
  
  for(i in 1:length(prop))
  {
    #MC1: makeConsensus classico con proportion da 0 a 1
    st1[i] <- system.time(ce <- clusterSingle(simData, subsample = TRUE, sequential = FALSE,
                                           mainClusterArgs = list(clusterFunction = "hierarchical01",
                                                                  clusterArgs = list(alpha = prop[i])),
                                           subsampleArgs = list(clusterFunction = "kmeans",
                                                                clusterArgs = list(k = 5),
                                                                samp.p = 0.7,
                                                                resamp.num = 100)))[[3]]
    nclust1[i]=length(table(ce@clusterMatrix))
    idx1[i]=adjustedRandIndex(trueCluster,ce@clusterMatrix)
    
  }
  #questa matrice è una per ognuno dei metodi, quindi 5 in totale
  mc=cbind(nclust1,st1,idx1)
  return(list(mc=mc))
}

ngruppi=10
range=c(-10.5,10.5)

sim=analisi(ngruppi,range)

nclust_mc=sim[[1]][,1]

st_mc=sim[[1]][,2]


idx_mc=sim[[1]][,3]


R=99

set.seed(123) 
for(i in 1:R)
{
  cat(i)
  sim=analisi(ngruppi=10,range=c(-5.5,5.5))
  nclust_mc=cbind(nclust_mc,sim[[1]][,1])
  
  st_mc=cbind(st_mc,sim[[1]][,2])
  
  idx_mc=cbind(idx_mc,sim[[1]][,3])
}

output=list("Numero di cluster"=nclust_mc,
            
            "System time"=st_mc,
            
            "Indice di Rand"=idx_mc)

save(output,file="Statistiche subsampling gerarchico.RData")

rm(list=ls())
load("Statistiche subsampling gerarchico.RData")

plot(rowMeans(output[[1]]),type="l")
plot(rowMeans(output[[2]]),type="l")
plot(rowMeans(output[[3]]),type="l")

prop=seq(0,0.9,by=0.1)
plot(rowMeans(output[[1]]),col=1,type="l",xlab="proportion",xaxt="n",main="Subsampling - Gerarchico",ylab="nclust")
axis(1, at=1:length(prop), labels=prop)
plot(rowMeans(output[[2]]),col=1,type="l",xlab="proportion",xaxt="n",main="Subsampling - Gerarchico",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
plot(rowMeans(output[[3]]),col=1,type="l",xlab="proportion",xaxt="n",main="Subsampling - Gerarchico",ylab="Indice corretto di Rand")
axis(1, at=1:length(prop), labels=prop)
