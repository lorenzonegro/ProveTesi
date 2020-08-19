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

ngruppi=3
range=c(-5.5,5.5)
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
                                        clusterArgs = list(alpha = 0.3)),
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
  nclust2=st2=idx2=rep(NA,length(prop))
  nclust3=st3=idx3=rep(NA,length(prop))
  nclust4=st4=idx4=rep(NA,length(prop))
  nclust5=st5=idx5=rep(NA,length(prop))
  
  
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
    
    #MC2
    st2[i] <- system.time(ce <- clusterSingle(simData, subsample = TRUE, sequential = FALSE,
                                              mainClusterArgs = list(clusterFunction = "snn",
                                                                     clusterArgs = list(alpha = prop[i],algorithm="walktrap")),
                                              subsampleArgs = list(clusterFunction = "kmeans",
                                                                   clusterArgs = list(k = 5),
                                                                   samp.p = 0.7,
                                                                   resamp.num = 100)))[[3]]
    nclust2[i]=length(table(ce@clusterMatrix))
    idx2[i]=adjustedRandIndex(trueCluster,ce@clusterMatrix)
    
    #MC3
    st3[i] <- system.time(ce <- clusterSingle(simData, subsample = TRUE, sequential = FALSE,
                                              mainClusterArgs = list(clusterFunction = "snn",
                                                                     clusterArgs = list(alpha = prop[i]),algorithm="louvain"),
                                              subsampleArgs = list(clusterFunction = "kmeans",
                                                                   clusterArgs = list(k = 5),
                                                                   samp.p = 0.7,
                                                                   resamp.num = 100)))[[3]]
    nclust3[i]=length(table(ce@clusterMatrix))
    idx3[i]=adjustedRandIndex(trueCluster,ce@clusterMatrix)
    
    #MC4
    st4[i] <- system.time(ce <- clusterSingle(simData, subsample = TRUE, sequential = FALSE,
                                              mainClusterArgs = list(clusterFunction = "snn",
                                                                     clusterArgs = list(alpha = prop[i]),algorithm="comp_weak"),
                                              subsampleArgs = list(clusterFunction = "kmeans",
                                                                   clusterArgs = list(k = 5),
                                                                   samp.p = 0.7,
                                                                   resamp.num = 100)))[[3]]
    nclust4[i]=length(table(ce@clusterMatrix))
    idx4[i]=adjustedRandIndex(trueCluster,ce@clusterMatrix)
    
    #MC5
    st5[i] <- system.time(ce <- clusterSingle(simData, subsample = TRUE, sequential = FALSE,
                                              mainClusterArgs = list(clusterFunction = "snn",
                                                                     clusterArgs = list(alpha = prop[i]),algorithm="comp_strong"),
                                              subsampleArgs = list(clusterFunction = "kmeans",
                                                                   clusterArgs = list(k = 5),
                                                                   samp.p = 0.7,
                                                                   resamp.num = 100)))[[3]]
    nclust5[i]=length(table(ce@clusterMatrix))
    idx5[i]=adjustedRandIndex(trueCluster,ce@clusterMatrix)
    
  }
  #questa matrice è una per ognuno dei metodi, quindi 5 in totale
  mc1=cbind(nclust1,st1,idx1)
  mc_walktrap=cbind(nclust2,st2,idx2)
  mc_louvain=cbind(nclust3,st3,idx3)
  mc_comp.weak=cbind(nclust4,st4,idx4)
  mc_comp.strong=cbind(nclust5,st5,idx5)
  return(list(mc=mc1,mc_walktrap=mc_walktrap,mc_louvain=mc_louvain,
              mc_comp.weak=mc_comp.weak,mc_comp.strong=mc_comp.strong))
}

ngruppi=3
range=c(-5.5,5.5)

sim=analisi(ngruppi,range)
prop=seq(0,0.9,by=0.1)

nclust_mc=sim[[1]][,1]
nclust_wa=sim[[2]][,1]
nclust_lo=sim[[3]][,1]
nclust_cw=sim[[4]][,1]
nclust_cs=sim[[5]][,1]

st_mc=sim[[1]][,2]
st_wa=sim[[2]][,2]
st_lo=sim[[3]][,2]
st_cw=sim[[4]][,2]
st_cs=sim[[5]][,2]

idx_mc=sim[[1]][,3]
idx_wa=sim[[2]][,3]
idx_lo=sim[[3]][,3]
idx_cw=sim[[4]][,3]
idx_cs=sim[[5]][,3]
R=99

#prima simulazione casuale
set.seed(123) #simulazioni 2-50
#set.seed(456) #simulazioni 51-100
for(i in 1:R)
{
  cat(i)
  sim=analisi(ngruppi=3,range=c(-5.5,5.5))
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


save(output,file="Statistiche subsampling gerarchico 3 vicini.RData")

rm(list=ls())

load("Statistiche subsampling gerarchico 3 vicini.RData")

nclust_mc=round(rowMeans(output$`Numero di cluster mC`))
nclust_wa=round(rowMeans(output$`Numero di cluster walktrap`))
nclust_lo=round(rowMeans(output$`Numero di cluster Louvain`))
nclust_cw=round(rowMeans(output$`Numero di cluster comp_weak`))
nclust_cs=round(rowMeans(output$`Numero di cluster comp_strong`))

prop=seq(0,0.9,by=0.1)

#confronto 
par(mfrow=c(1,3))
plot(nclust_mc,col=1,type="l",xlab="proportion",xaxt="n",main="Gerarchico",ylab="Numero di cluster")
axis(1, at=1:length(prop), labels=prop)
plot(nclust_wa,col=2,type="l",xlab="proportion",xaxt="n",main="Louvain & walktrap",ylab="Numero di cluster")
axis(1, at=1:length(prop), labels=prop)
points(nclust_lo,col=3,type="l")
legend("topleft",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(nclust_cw,col=4,type="l",xlab="proportion",xaxt="n",main="Components",ylab="Numero di cluster")
axis(1, at=1:length(prop), labels=prop)
points(nclust_cs,col=6,type="l")
legend("topleft",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)

#indice
idx_mc=rowMeans(output$`Indice mc`)
idx_wa=rowMeans(output$`Indice walktrap`)
idx_lo=rowMeans(output$`Indice Louvain`)
idx_cw=rowMeans(output$`Indice comp_weak`)
idx_cs=rowMeans(output$`Indice comp_strong`)

par(mfrow=c(1,3))
plot(idx_mc,col=1,type="l",xlab="proprtion",xaxt="n",main="Gerarchico",ylab="Adjusted Rand Index")
axis(1, at=1:length(prop), labels=prop)
plot(idx_wa,col=2,type="l",xlab="proportion",xaxt="n",main="Walktrap & Louvain",ylab="Adjusted Rand Index")
axis(1, at=1:length(prop), labels=prop)
points(idx_lo,col=3,type="l")
legend("topright",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.5)
plot(idx_cw,col=4,type="l",xlab="proportion",xaxt="n",main="Components",ylab="Adjusted Rand Index")
axis(1, at=1:length(prop), labels=prop)
points(idx_cs,col=6,type="l")
legend("topright",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.5)

#system time
st_mc=rowMeans(output$`System time mC`)
st_wa=rowMeans(output$`System time walktrap`)
st_lo=rowMeans(output$`System time Louvain`)
st_cw=rowMeans(output$`System time comp_weak`)
st_cs=rowMeans(output$`System time comp_strong`)

par(mfrow=c(1,3))
plot(st_mc,col=1,type="l",xlab="proprtion",xaxt="n",main="Gerarchico",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
plot(st_wa,col=2,type="l",xlab="proportion",xaxt="n",main="Walktrap & Louvain",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
points(st_lo,col=3,type="l")
legend("topright",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.5)
plot(st_cw,col=4,type="l",xlab="proportion",xaxt="n",main="Components",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
points(st_cs,col=6,type="l")
legend("topright",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.5)

#Confronto best performance
best_mc=mean(apply(output$`Indice mc`,2,function(x) max(x)))
idx_best_mc=apply(output$`Indice mc`,2,function(x) which.max(x))
st_best_mc=mean(st_mc[idx_best_mc])

best_wa=mean(apply(output$`Indice walktrap`,2,function(x) max(x)))
idx_best_wa=apply(output$`Indice walktrap`,2,function(x) which.max(x))
st_best_wa=mean(st_wa[idx_best_wa])

best_lo=mean(apply(output$`Indice Louvain`,2,function(x) max(x)))
idx_best_lo=apply(output$`Indice Louvain`,2,function(x) which.max(x))
st_best_lo=mean(st_lo[idx_best_lo])

best_cw=mean(apply(output$`Indice comp_weak`,2,function(x) max(x)))
idx_best_cw=apply(output$`Indice comp_weak`,2,function(x) which.max(x))
st_best_cw=mean(st_cw[idx_best_cw])

best_cs=mean(apply(output$`Indice comp_strong`,2,function(x) max(x)))
idx_best_cs=apply(output$`Indice comp_strong`,2,function(x) which.max(x))
st_best_cs=mean(st_cs[idx_best_cs])

# Libraries
library(ggplot2)

# Create data
data <- data.frame(
  x=c("Gerarchico","Walktrap","Louvain","Comp_Weak","Comp_Strong"),
  y=c(best_mc,best_wa,best_lo,best_cw,best_cs)
)

# Horizontal version
ggplot(data, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  labs(x="Algoritmh",y="Adjusted Rand Index")+
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# Create data
data <- data.frame(
  x=c("Gerarchcio","Walktrap","Louvain","Comp_Weak","Comp_Strong"),
  y=c(st_best_mc,st_best_wa,st_best_lo,st_best_cw,st_best_cs)
)

# Horizontal version
ggplot(data, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  labs(x="Algoritmh",y="System time")+
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


table(prop[idx_best_mc])
table(prop[idx_best_wa])
table(prop[idx_best_lo])
table(prop[idx_best_cw])
table(prop[idx_best_cs])
