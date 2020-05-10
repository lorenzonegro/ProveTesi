####SimData Esteso 10 gruppi non omogenei####
nvar<-51 #multiple of 3
set.seed(1000)
num_grp <- rep(0,10)
for(i in 1:10)
{
  num_grp[i]=round(runif(1,10,600))
}
num_grp=num_grp-round((sum(num_grp)-3000)/10)
sum(num_grp)
caso=round(runif(1,0.5,10.5))
num_grp[caso]=num_grp[caso]+1
sum(num_grp)

set.seed(123)
medie=runif(10,-5.5,5.5)
x<-cbind(matrix(rnorm(num_grp[1]*nvar,mean=medie[1]),nrow=nvar),
         matrix(rnorm(num_grp[2]*nvar,mean=medie[2]),nrow=nvar),
         matrix(rnorm(num_grp[3]*nvar,mean=medie[3]),nrow=nvar),
         matrix(rnorm(num_grp[4]*nvar,mean=medie[4]),nrow=nvar),
         matrix(rnorm(num_grp[5]*nvar,mean=medie[5]),nrow=nvar),
         matrix(rnorm(num_grp[6]*nvar,mean=medie[6]),nrow=nvar),
         matrix(rnorm(num_grp[7]*nvar,mean=medie[7]),nrow=nvar),
         matrix(rnorm(num_grp[8]*nvar,mean=medie[8]),nrow=nvar),
         matrix(rnorm(num_grp[9]*nvar,mean=medie[9]),nrow=nvar),
         matrix(rnorm(num_grp[10]*nvar,mean=medie[10]),nrow=nvar))
#make some of them flipped effects (better for testing if both sig under/over
#expressed variables)
prob=num_grp/sum(num_grp)
geneGroup<-sample(1:10,size=50,prob=prob,replace=T)
gpIndex<-list(1:num_grp[1],(num_grp[1]+1):(num_grp[2]),
              (num_grp[2]+1):(num_grp[3]),(num_grp[3]+1):(num_grp[4]),
              (num_grp[4]+1):(num_grp[5]),(num_grp[5]+1):(num_grp[6]),
              (num_grp[6]+1):(num_grp[7]),(num_grp[7]+1):(num_grp[8]),
              (num_grp[8]+1):(num_grp[9]),(num_grp[9]+1):(num_grp[10]))


#add in differences in variable means
smp<-sample(1:nrow(x),10)
x[smp,]<-x[smp,]+10

#make different signal y
set.seed(456)
medie=runif(10,-1,1)
y<-cbind(matrix(rnorm(num_grp[1]*nvar,mean=medie[1]),nrow=nvar),
         matrix(rnorm(num_grp[2]*nvar,mean=medie[2]),nrow=nvar),
         matrix(rnorm(num_grp[3]*nvar,mean=medie[3]),nrow=nvar),
         matrix(rnorm(num_grp[4]*nvar,mean=medie[4]),nrow=nvar),
         matrix(rnorm(num_grp[5]*nvar,mean=medie[5]),nrow=nvar),
         matrix(rnorm(num_grp[6]*nvar,mean=medie[6]),nrow=nvar),
         matrix(rnorm(num_grp[7]*nvar,mean=medie[7]),nrow=nvar),
         matrix(rnorm(num_grp[8]*nvar,mean=medie[8]),nrow=nvar),
         matrix(rnorm(num_grp[9]*nvar,mean=medie[9]),nrow=nvar),
         matrix(rnorm(num_grp[10]*nvar,mean=medie[10]),nrow=nvar))
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