rm(list=ls())
library(devtools)
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/BiocNeighbors")
load_all()
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/clusterExperiment")
load_all()
library(scran)
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/ProveTesi")

####SimData Esteso 10 gruppi####
nvar<-51 #multiple of 3
n<-300
set.seed(123)
medie=runif(10,-5.5,5.5)
x<-cbind(matrix(rnorm(n*nvar,mean=medie[1]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[2]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[3]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[4]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[5]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[6]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[7]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[8]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[9]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[10]),nrow=nvar))
#make some of them flipped effects (better for testing if both sig under/over
#expressed variables)
geneGroup<-sample(rep(1:10,each=floor(nvar/10)))
gpIndex<-list(1:n,(n+1):(n*2),
              (n*2+1):(n*3),(n*3+1):(n*4),
              (n*4+1):(n*5),(n*5+1):(n*6),
              (n*6+1):(n*7),(n*7+1):(n*8),
              (n*8+1):(n*9),(n*9+1):(n*10))
x[geneGroup==1,]<-x[geneGroup==1,unlist(gpIndex[c(10,1,2,3,4,5,6,7,8,9)])]
x[geneGroup==2,]<-x[geneGroup==2,unlist(gpIndex[c(9,10,1,2,3,4,5,6,7,8)])]
x[geneGroup==3,]<-x[geneGroup==3,unlist(gpIndex[c(8,9,10,1,2,3,4,5,6,7)])]
x[geneGroup==4,]<-x[geneGroup==4,unlist(gpIndex[c(7,8,9,10,1,2,3,4,5,6)])]
x[geneGroup==5,]<-x[geneGroup==5,unlist(gpIndex[c(6,7,8,9,10,1,2,3,4,5)])]
x[geneGroup==6,]<-x[geneGroup==6,unlist(gpIndex[c(5,6,7,8,9,10,1,2,3,4)])]
x[geneGroup==7,]<-x[geneGroup==7,unlist(gpIndex[c(4,5,6,7,8,9,10,1,2,3)])]
x[geneGroup==8,]<-x[geneGroup==8,unlist(gpIndex[c(3,4,5,6,7,8,9,10,1,2)])]
x[geneGroup==9,]<-x[geneGroup==9,unlist(gpIndex[c(2,3,4,5,6,7,8,9,10,1)])]

#add in differences in variable means
smp<-sample(1:nrow(x),10)
x[smp,]<-x[smp,]+10

#make different signal y
set.seed(456)
medie=runif(10,-1,1)
y<-cbind(matrix(rnorm(n*nvar,mean=medie[1]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[2]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[3]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[4]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[5]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[6]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[7]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[8]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[9]),nrow=nvar),
         matrix(rnorm(n*nvar,mean=medie[10]),nrow=nvar))
y<-y[,sample(1:ncol(y))]+ matrix(rnorm(10*n*nvar,sd=3),nrow=nvar)

#add together the two signals
simData<-x+y

#add pure noise variables
simData<-rbind(simData,matrix(rnorm(10*n*nvar,mean=10),nrow=nvar),
               matrix(rnorm(10*n*nvar,mean=5),nrow=nvar))
#make count data
countMean<-exp(simData/2)
simCount<-matrix(rpois(n=length(as.vector(countMean)), lambda
                       =as.vector(countMean)+.1),nrow=nrow(countMean),ncol=ncol(countMean))
#labels for the truth
trueCluster<-rep(c(1:10),each=n)

save(simData,file="simData_10.RData")
save(simCount,file="simCount_10.RData")
save(trueCluster,file="trueCluster_10.RData")

####Analisi####
cl=clusterMany(simData,k=2:20,clusterFunction="kmeans")@clusterMatrix #k numero cluster

mc1=makeConsensus(cl,proportion=0.9)
mc2=makeConsensus2(cl,k=10,algorithm="cluster_walktrap") #aggiungere step a makeconsensus2
mc3=makeConsensus2(cl,k=10,alghorithm="components_strong")

#The weakly connected components are found by a simple breadth-first search.
#The strongly connected components are implemented by two consecutive depth-first searches.

table(mc1$clustering)
table(mc2$clustering)
table(mc3$clustering)
table(trueCluster)

#aggiungi system.time
system.time(makeConsensus(cl,proportion=0.7))
system.time(makeConsensus2(cl,k=10,algorithm="cluster_louvain"))

Rprof(filename ="ram_consensus.txt", append = FALSE, memory.profiling = TRUE)
mc1=makeConsensus(cl,proportion=0.7)
Rprof(NULL)

profile <- summaryRprof(filename = "ram_consensus.txt", chunksize = -1L, 
                        memory = "tseries", diff = FALSE)
max_mem_cons1 <- max(rowSums(profile[,1:3]))*0.00000095367432
max_mem_cons1

Rprof(filename ="ram_consensus2.txt", append = FALSE, memory.profiling = TRUE)
mc2=makeConsensus2(cl,k=10,algorithm="cluster_walktrap")
Rprof(NULL)

profile2 <- summaryRprof(filename = "ram_consensus2.txt", chunksize = -1L, 
                         memory = "tseries", diff = FALSE)
max_mem_cons2 <- max(rowSums(profile2[,1:3]))*0.00000095367432
max_mem_cons2

#####Merge Cluster####
cl=clusterMany(simData,k=2:10,clusterFunction="kmeans")
ce<-makeConsensus2(cl, k=10, algorithm="cluster_walktrap")
ce<-makeDendrogram(ce,reduceMethod="mad") #pca al posto di mad
ce<-mergeClusters(ce,mergeMethod="adjP",DEMethod="limma",cutoff=0.3,plot=FALSE)
ce
#come gestire?
#e come gestire makeConsensus2

#clusterFunction ?
#guardare indice
#mc2 abbiamo due algoritmi mc1 ne abbiamo uno vogliamo confrontarli

#lascio stare per ora merge cluster
#per ogni metodo 3 grafici
#asse x cambio parametro (k mc2 proprtion mc1)
#asse y, 4 grafici, system.time (memory fallo successivamente) justrandindex, numero cluster
#10 vaori per proportion 10 valori per k
