rm(list=ls())
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/ProveTesi")

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


medie=runif(3,-5.5,5.5)
prop=c(0.3,0.5,0.2)

simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=150,
                           k=3, x_mus = c(0,5,5), 
                           x_sds = c(1,0.1,1), 
                           y_mus = c(5,5,0), 
                           y_sds = c(1,0.1,1), 
                           prop1 = c(0.3,0.5,0.2))
str(simulazioni)
trueCluster=simulazioni$true_cluster_id
table(trueCluster)
simData=as.matrix(t(simulazioni$obs_data))
dim(simData)
save(simData,file="simData_3_casuali.RData")
save(trueCluster,file="trueCluster_3_casuali.RData")

medie=runif(10,-5.5,5.5)
prop=c(runif(10,0,1))
prop=prop/sum(prop)
prop

simulazioni=simulate_gauss_mix(n_cells=3000, n_genes=150,
                               k=10, x_mus = medie, 
                               x_sds = c(runif(10,0.1,2)), 
                               y_mus = c(5,5,0,5,5,0,5,5,0,5), 
                               y_sds = c(rep(1,10)), 
                               prop1 = prop)
str(simulazioni)
trueCluster=simulazioni$true_cluster_id
table(trueCluster)
simData=as.matrix(t(simulazioni$obs_data))
dim(simData)
save(simData,file="simData_10_casuali.RData")
save(trueCluster,file="trueCluster_3_casuali.RData")
