sim=analisi(ngruppi=10,range=c(-5.5,5.5))
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
