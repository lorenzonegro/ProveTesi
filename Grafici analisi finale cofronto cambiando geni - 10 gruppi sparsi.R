rm(list=ls())

load("Statistiche cambio numero geni.RData")

numero_geni=c(50,150,300,500,1000)

#indice
idx_mc=c(output$`Indice mc`)
idx_lo=c(output$`Indice Louvain`)

plot(idx_mc,col=1,type="l",xlab="Numero di geni",xaxt="n",main="mC vs mC3(Louvain)",ylab="Adjusted Rand Index",ylim=c(min(idx_mc,idx_lo),1))
axis(1, at=1:length(numero_geni), labels=numero_geni)
points(idx_lo,col=3,type="l")
legend("bottomleft",legend=c("mC","mC3"),col=c("black","green"),lty=1,cex=0.5)

#system time
st_mc=c(output$`System time mC`)
st_lo=c(output$`System time Louvain`)

plot(st_mc,col=1,type="l",xlab="Numero di geni",xaxt="n",main="mC vs mC3(Louvain)",ylab="System time",ylim=c(min(st_lo),max(st_mc)))
axis(1, at=1:length(numero_geni), labels=numero_geni)
points(st_lo,col=3,type="l")
legend("bottomleft",legend=c("mC","mC3"),col=c("black","green"),lty=1,cex=0.5)

#system time
mem_mc=c(output$`Memoria mC`)
mem_lo=c(output$`Memoria mC3(Louvain)`)

plot(mem_mc,col=1,type="l",xlab="Numero di geni",xaxt="n",main="mC vs mC3(Louvain)",ylab="Consumo memoria",ylim=c(min(mem_lo),max(mem_mc)))
axis(1, at=1:length(numero_geni), labels=numero_geni)
points(mem_lo,col=3,type="l")
legend("topleft",legend=c("mC","mC3"),col=c("black","green"),lty=1,cex=0.8)
