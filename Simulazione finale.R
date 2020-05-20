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

prop=seq(0,0.9,by=0.1)


analisi <- function()
{
  nclust1=st1=idx1=rep(NA,length(prop))
  nclust2=st2=idx2=rep(NA,length(prop))
  nclust3=st3=idx3=rep(NA,length(prop))
  nclust4=st4=idx4=rep(NA,length(prop))
  nclust5=st5=idx5=rep(NA,length(prop))
  
  cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix
  adj_CM=apply(cl,2,function(x) adjustedRandIndex(trueCluster,x))
  
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

    #MC3: makeConsensus nuovo prop diversi algoritmo walktrap
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
  mc_louvain=cbind(nclust1,st1,idx1)
  mc_walktrap=cbind(nclust1,st1,idx1)
  mc_comp.weak=cbind(nclust1,st1,idx1)
  mc_comp.strong=cbind(nclust1,st1,idx1)
  return(list(mc1,mc_louvain,mc_walktrap,mc_comp.weak,mc_comp.strong))
}