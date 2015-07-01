library(ergm)
library(Rmpi)
library(snow)
library(influenceR)

data(flo)
fe <- which(flo>0, arr.ind=T)
fg <- igraph::graph(fe, directed=F)

cl <- makeMPIcluster(2)

bridging(fg, MPI=T, cluster=cl)

stopCluster(cl)
mpi.exit()

