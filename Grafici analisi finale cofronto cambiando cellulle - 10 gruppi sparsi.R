rm(list=ls())

load("Statistiche cambio numero cellule.RData")

numero_cellule=c(1000,3000,5000,7000,10000)

#indice
idx_mc=c(output$`Indice mc`)
idx_lo=c(output$`Indice Louvain`)

plot(idx_mc,col=1,type="l",xlab="Numero di Cellule",xaxt="n",main="mC vs mC3(Louvain)",ylab="Adjusted Rand Index",ylim=c(min(idx_mc,idx_lo),1))
axis(1, at=1:length(numero_cellule), labels=numero_cellule)
points(idx_lo,col=3,type="l")
legend("bottomleft",legend=c("mC","mC3"),col=c("black","green"),lty=1,cex=0.8)

#system time
st_mc=c(output$`System time mC`)
st_lo=c(output$`System time Louvain`)

plot(st_mc,col=1,type="l",xlab="Numero di Cellule",xaxt="n",main="mC vs mC3(Louvain)",ylab="System time")
axis(1, at=1:length(numero_cellule), labels=numero_cellule)
points(st_lo,col=3,type="l")
legend("topleft",legend=c("mC","mC3"),col=c("black","green"),lty=1,cex=0.8)

#system time
mem_mc=c(output$`Memoria mC`)
mem_lo=c(output$`Memoria mC3(Louvain)`)

plot(mem_mc,col=1,type="l",xlab="Numero di Cellule",xaxt="n",main="mC vs mC3(Louvain)",ylab="Consumo memoria")
axis(1, at=1:length(numero_cellule), labels=numero_cellule)
points(mem_lo,col=3,type="l")
legend("topleft",legend=c("mC","mC3"),col=c("black","green"),lty=1,cex=0.8)
