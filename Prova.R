nobs <- 10000
ndim <- 20
data <- matrix(runif(nobs*ndim), ncol=ndim)
library(BiocNeighbors)
fout <- findKNN(data, k=10, BNPARAM=KmknnParam())
head(fout$index)

#Esempio 1
#10 osservazioni 2 dimensioni
data1 <- matrix(c(1,1,10,10,1,2,20,20,2,1,20,30,30,30,300,100,200,100,0.3,0.3),ncol=2,byrow=T)
data1
fout1 <- findKNN(data1, k=2, BNPARAM=KmknnParam()) #k=2 2 vicini più vicini
head(fout1$index) #Mi dice quali sono i vicni più vicini
head(fout1$distance) #mi dice la sitanza da questi vicini
#Esempio
data1[1,] #ha come vicin più vicini
head(fout1$index) #10 e 3
data1[10,]
data1[3,]
#la sua distanza da 10 e 3 è
head(fout1$distance)[1,]

#al posto di findKNN posso usare queryKNN per confrontare due dataset
query1 <- matrix(c(10,10,1,1,10,20,2,2,20,10,2,3,3,3,0.3,1,200,10,0.3,0.3),ncol=2,byrow=T)
query1
#Confronto data1 e query1
qout <- queryKNN(data1, query1, k=2, BNPARAM=KmknnParam(distance="Euclidean")) #distance può essere manhattan o euclidean
head(qout$index) #I punti più vicni a primo punto di query sono il secondo e il terzo in data
query1[1,]
data1[2,]
data1[3,]

#Esempio 2 x Hamming Distance
#9 osservazioni 3 cluster e 3 clustering costruiamo una matrice 6x3
#k=3 3 vicini più vicini
data2 <- matrix(c(1,2,3,1,1,1,2,2,2,1,2,3,3,3,3,1,2,1,3,3,3,1,1,2,2,2,3),ncol=3,byrow=T)
data2
fout2 <- findKNN(data2,k=3,BNPARAM=KmknnParam())
fout2$index

#Come implementare la distanza di Hamming con R?
#Esiste una funzione Hamming Distance ma funziona in modo booleano
#Quindi se due elementi sono zero (o entrambi diversi da 0) mi da true
#Altrimenti se uno è 0 e l'altro no mi da false
#Io voglio implementare una funzione che confronti i singoli elementi
hamming.distance<-function (x, y) 
{
  x=paste(x,collapse="")
  y=paste(y,collapse="")
  
  if(nchar(x)!=nchar(y))
  {
    stop("x and y must have the same length")
  }
  z <- adist(x,y)[1,1] #Misura la distanza tra stringhe
  return(z)
}

x=c(3,4,0,1)
y=c(1,1,1,1)
hamming.distance(x,y)

x="Ciao"
y="Ciad"
hamming.distance(x,y)

x=c("C","i","a","o")
y=c("C","i","a","d")
hamming.distance(x,y)

y="Ciao"
x="Coia"
hamming.distance(x,y) #perché viene 2? Dovrebbe essere 3

#Sembra funzionare, tranne in quel caso. Una versione meno "smart" invece sembra funzionare sempre
hamming.distance<-function (x, y) 
{
  x=paste(x,collapse="")
  y=paste(y,collapse="")
  if(nchar(x)!=nchar(y))
  {
    stop("x and y must have the same length")
  }
  l=nchar(x)
  z <- apply(do.call(rbind, strsplit(c(x, y), "")), 2, 
             function(x){length(unique(x[!x %in% "_"])) == 1})
  
  z=l-sum(z)
  return(z)
}

x=c(3,4,0,1)
y=c(1,1,1,1)
hamming.distance(x,y)

x="Ciao"
y="Ciad"
hamming.distance(x,y)

x=c("C","i","a","o")
y=c("C","i","a","d")
hamming.distance(x,y)

y="Ciao"
x="Coia"
hamming.distance(x,y) #Ora viene 3

