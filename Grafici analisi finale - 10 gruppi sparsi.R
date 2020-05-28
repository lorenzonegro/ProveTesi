rm(list=ls())

load("Statistiche 10 gruppi sparsi.RData")

nclust_mc=round(rowMeans(output$`Numero di cluster mC`))
nclust_wa=round(rowMeans(output$`Numero di cluster walktrap`))
nclust_lo=round(rowMeans(output$`Numero di cluster Louvain`))
nclust_cw=round(rowMeans(output$`Numero di cluster walktrap`))
nclust_cs=round(rowMeans(output$`Numero di cluster comp_strong`))

prop=seq(0,0.9,by=0.1)

#confronto 
par(mfrow=c(1,3))
plot(nclust_mc,col=1,type="l",xlab="proportion",xaxt="n",main="makeConsensus",ylab="Numero di cluster")
axis(1, at=1:length(prop), labels=prop)
plot(nclust_wa,col=2,type="l",xlab="proportion",xaxt="n",main="mc3 - Louvain & walktrap",ylab="Numero di cluster")
axis(1, at=1:length(prop), labels=prop)
points(nclust_lo,col=3,type="l")
legend("topleft",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(nclust_cw,col=4,type="l",xlab="proportion",xaxt="n",main="mc3 - components",ylab="Numero di cluster")
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
plot(idx_mc,col=1,type="l",xlab="proprtion",xaxt="n",main="makeConsensus",ylab="Adjusted Rand Index")
axis(1, at=1:length(prop), labels=prop)
plot(idx_wa,col=2,type="l",xlab="proportion",xaxt="n",main="mc3 - Walktrap & Louvain",ylab="Adjusted Rand Index")
axis(1, at=1:length(prop), labels=prop)
points(idx_lo,col=3,type="l")
legend("bottomleft",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(idx_cw,col=4,type="l",xlab="proportion",xaxt="n",main="mc3 - components",ylab="Adjusted Rand Index")
axis(1, at=1:length(prop), labels=prop)
points(idx_cs,col=6,type="l")
legend("bottomleft",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)

#system time
st_mc=rowMeans(output$`System time mC`)
st_wa=rowMeans(output$`System time walktrap`)
st_lo=rowMeans(output$`System time Louvain`)
st_cw=rowMeans(output$`System time comp_weak`)
st_cs=rowMeans(output$`System time comp_strong`)

par(mfrow=c(1,3))
plot(st_mc,col=1,type="l",xlab="proprtion",xaxt="n",main="makeConsensus",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
plot(st_wa,col=2,type="l",xlab="proportion",xaxt="n",main="mC3 - Walktrap & Louvain",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
points(st_lo,col=3,type="l")
legend("topright",legend=c("walktrap","louvain"),col=c("red","green"),lty=1,cex=0.8)
plot(st_cw,col=4,type="l",xlab="proportion",xaxt="n",main="mC3 - components",ylab="System time")
axis(1, at=1:length(prop), labels=prop)
points(st_cs,col=6,type="l")
legend("topright",legend=c("weak","strong"),col=c("blue","violet"),lty=1,cex=0.8)

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
  x=c("mC","Walktrap","Louvain","Comp_Weak","Comp_Strong"),
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
  x=c("mC","Walktrap","Louvain","Comp_Weak","Comp_Strong"),
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
