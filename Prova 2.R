library(devtools)
load_all()
data=matrix(c(1,1,1,
              1,1,2,
              3,1,2,
              2,2,2,
              1,1,1),5,3,byrow=T)

clusterMatrix(data)

fout2 <- findKNN(data,k=2,BNPARAM=VptreeParam(distance="Hamming"))
fout2$index
data
