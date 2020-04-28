#Prova 6 
#prova di utilizzo di Merge cluster

#Riprendo prova 4
#Aggiungo Mergecluster
#Valutare anche system time e uso della RAM
#Rprof x utilizzo ram (vedi skype)
#simulare da normale trivariata in stile ?simData 3000 dati al posto di 300
#adjusted rand index nel pachetto mclust anche se il numero di cluster non corrisponde

####Preliminare####
rm(list=ls())
library(devtools)
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/BiocNeighbors")
load_all()
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/clusterExperiment")
load_all()
library(scran)
setwd("C:/Users/user1/Desktop/UNIVERSITA'/TESI/ProveTesi")

####SimData Esteso####
nvar<-51 #multiple of 3
n<-1000
x<-cbind(matrix(rnorm(n*nvar,mean=5),nrow=nvar),
         matrix(rnorm(n*nvar,mean=-5),nrow=nvar),
         matrix(rnorm(n*nvar,mean=0),nrow=nvar))
#make some of them flipped effects (better for testing if both sig under/over
#expressed variables)
geneGroup<-sample(rep(1:3,each=floor(nvar/3)))
gpIndex<-list(1:n,(n+1):(n*2),(2*n+1):(n*3))
x[geneGroup==1,]<-x[geneGroup==1,unlist(gpIndex[c(3,1,2)])]
x[geneGroup==2,]<-x[geneGroup==2,unlist(gpIndex[c(2,3,1)])]

#add in differences in variable means
smp<-sample(1:nrow(x),10)
x[smp,]<-x[smp,]+10

#make different signal y
y<-cbind(matrix(rnorm(n*nvar,mean=1),nrow=nvar),
         matrix(rnorm(n*nvar,mean=-1),nrow=nvar),
         matrix(rnorm(n*nvar,mean=0),nrow=nvar))
y<-y[,sample(1:ncol(y))]+ matrix(rnorm(3*n*nvar,sd=3),nrow=nvar)

#add together the two signals
simData<-x+y

#add pure noise variables
simData<-rbind(simData,matrix(rnorm(3*n*nvar,mean=10),nrow=nvar),
               matrix(rnorm(3*n*nvar,mean=5),nrow=nvar))
#make count data
countMean<-exp(simData/2)
simCount<-matrix(rpois(n=length(as.vector(countMean)), lambda
                       =as.vector(countMean)+.1),nrow=nrow(countMean),ncol=ncol(countMean))
#labels for the truth
trueCluster<-rep(c(1:3),each=n)

####Analisi####
cl=clusterMany(simData,k=2:10,clusterFunction="kmeans")@clusterMatrix #k numero cluster

mc1=makeConsensus(cl,proportion=0.7)
mc2=makeConsensus2(cl,k=10,algorithm="cluster_walktrap") #aggiungere step a makeconsensus2

table(mc1$clustering)
table(mc2$clustering)
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
