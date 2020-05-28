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

analisi <- function(ngruppi,range,n_cells=3000,n_genes=150)
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
  
  prop=seq(0,0.9,by=0.1)
  nclust1=st1=idx1=rep(NA,length(prop))
  nclust2=st2=idx2=rep(NA,length(prop))
  nclust3=st3=idx3=rep(NA,length(prop))
  nclust4=st4=idx4=rep(NA,length(prop))
  nclust5=st5=idx5=rep(NA,length(prop))
  
  cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix
  
  for(i in 1:length(prop))
  {
    #MC1: makeConsensus classico con proportion da 0 a 1
    st1[i] <- system.time(mc1 <- makeConsensus(cl,proportion=prop[i]))[[3]]
    nclust1[i]=length(table(mc1$clustering))
    idx1[i]=adjustedRandIndex(trueCluster,mc1$clustering)
    
    #MC2: makeConsensus nuovo prop diversi algoritmo walktrap
    st2[i] <- system.time(mc2 <- makeConsensus3(cl, proportion=prop[i], algorithm="cluster_walktrap"))[[3]]
    nclust2[i]=length(table(mc2$clustering))
    idx2[i]=adjustedRandIndex(trueCluster,mc2$clustering)

    #MC3: makeConsensus nuovo prop diversi algoritmo louvain
    st3[i] <- system.time(mc3 <- makeConsensus3(cl, proportion=prop[i], algorithm="cluster_louvain"))[[3]]
    nclust3[i]=length(table(mc3$clustering))
    idx3[i]=adjustedRandIndex(trueCluster,mc3$clustering)
  
    st4[i] <- system.time(mc4 <- makeConsensus3(cl, proportion=prop[i], algorithm="components_weak"))[[3]]
    nclust4[i]=length(table(mc4$clustering))
    idx4[i]=adjustedRandIndex(trueCluster,mc4$clustering)

    st5[i] <- system.time(mc5 <- makeConsensus3(cl, proportion=prop[i], algorithm="components_strong"))[[3]]
    nclust5[i]=length(table(mc5$clustering))
    idx5[i]=adjustedRandIndex(trueCluster,mc5$clustering)
  }
  #restituisco in output una matrice con numero di righe pari a proportion (10)
  #3 colonne, una per system time, una per nclust, una per l'indice
  #questa matrice è una per ognuno dei metodi, quindi 5 in totale
  mc1=cbind(nclust1,st1,idx1)
  mc_walktrap=cbind(nclust2,st2,idx2)
  mc_louvain=cbind(nclust3,st3,idx3)
  mc_comp.weak=cbind(nclust4,st4,idx4)
  mc_comp.strong=cbind(nclust5,st5,idx5)
  return(list(mc=mc1,mc_walktrap=mc_walktrap,mc_louvain=mc_louvain,
              mc_comp.weak=mc_comp.weak,mc_comp.strong=mc_comp.strong))
}

load("Statistiche 10 gruppi.RData")

sim=output

prop=seq(0,0.9,by=0.1)

nclust_mc=sim[[1]]
nclust_wa=sim[[2]]
nclust_lo=sim[[3]]
nclust_cw=sim[[4]]
nclust_cs=sim[[5]]

st_mc=sim[[6]]
st_wa=sim[[7]]
st_lo=sim[[8]]
st_cw=sim[[9]]
st_cs=sim[[10]]

idx_mc=sim[[11]]
idx_wa=sim[[12]]
idx_lo=sim[[13]]
idx_cw=sim[[14]]
idx_cs=sim[[15]]

R=50

#prima simulazione casuale
#set.seed(123) #simulazioni 2-50
set.seed(456) #simulazioni 51-100
for(i in 1:R)
{
  cat(i)
  sim=analisi(ngruppi=10,range=c(-5.5,5.5))
  nclust_mc=cbind(nclust_mc,sim[[1]][,1])
  nclust_wa=cbind(nclust_wa,sim[[2]][,1])
  nclust_lo=cbind(nclust_lo,sim[[3]][,1])
  nclust_cw=cbind(nclust_cw,sim[[4]][,1])
  nclust_cs=cbind(nclust_cs,sim[[5]][,1])

  st_mc=cbind(st_mc,sim[[1]][,2])
  st_wa=cbind(st_wa,sim[[2]][,2])
  st_lo=cbind(st_lo,sim[[3]][,2])
  st_cw=cbind(st_cw,sim[[4]][,2])
  st_cs=cbind(st_cs,sim[[5]][,2])

  idx_mc=cbind(idx_mc,sim[[1]][,3])
  idx_wa=cbind(idx_wa,sim[[2]][,3])
  idx_lo=cbind(idx_lo,sim[[3]][,3])
  idx_cw=cbind(idx_cw,sim[[4]][,3])
  idx_cs=cbind(idx_cs,sim[[5]][,3])
}

output=list("Numero di cluster mC"=nclust_mc,
            "Numero di cluster walktrap"=nclust_wa,
            "Numero di cluster Louvain"=nclust_lo,
            "Numero di cluster comp_weak"=nclust_cw,
            "Numero di cluster comp_strong"=nclust_cs,
            "System time mC"=st_mc,
            "System time walktrap"=st_wa,
            "System time Louvain"=st_lo,
            "System time comp_weak"=st_cw,
            "System time comp_strong"=st_cs,
            "Indice mc"=idx_mc,
            "Indice walktrap"=idx_wa,
            "Indice Louvain"=idx_lo,
            "Indice comp_weak"=idx_cw,
            "Indice comp_strong"=idx_cs)

save(output,file="Statistiche 10 gruppi.RData")

plot(rowMeans(output[[11]]),type="l")
plot(rowMeans(output[[12]]),type="l")
plot(rowMeans(output[[13]]),type="l")
plot(rowMeans(output[[14]]),type="l")
plot(rowMeans(output[[15]]),type="l")

