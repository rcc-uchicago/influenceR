# igraph package has arpack functions
# evcent - eigenvector centrality
# betweenness - vertex betweenness

dimacs.to.graph <- function(fname) {
    x <- read.table(fname, skip=1)
    el <- as.matrix(x[2:3])
    
    graph.edgelist(el, directed=F)
}

# requires network to be loaded
igraph.to.network <- function(g) {
  el <- igraph::get.edgelist(g)
  net <- network::as.network(el, directed=F) # directed?
  net
}

eigencentrality <- function(g) {
  evcent(g, scale=F)$vector
}

# 1/2 the values of our betweenness code, which is because this is UNDIRECTED for real
betweenness <- function(g, snap=T) {
  
  if (!snap)
    return(igraph::betweenness(g, directed=F))
  
  el <- get.edgelist(g)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el))
  m <- as.integer(length(el)/2) # TODO: for directed too?
  
  .Call("snap_betweenness_R", el_i, n, m, PACKAGE="influenceR")
}

# g : graph, k : size of set, prob : probability we will accept a lower state,
# tol : acceptance tolerance, maxsec : total computation budget,
# roundsec : seconds for a round (in parallel version)
keyplayer <- function(g, k, prob=0.0, tol=0.0001, maxsec=600, roundsec=30) {
  el <- get.edgelist(g)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el))
  m <- as.integer(length(el)/2)

  s <- .Call("snap_keyplayer_R", el_i, n, m, as.integer(k), prob, tol, as.integer(maxsec), as.integer(roundsec), PACKAGE="influenceR")
  
  which(s>0)
}


bridging <- function(g, MPI=F) {
  
  el <- get.edgelist(g)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el_i))
  m <- as.integer(length(el_i)/2)
  
  if(!MPI) {
    x <- .Call("snap_bridging_R", el_i, n, m, as.integer(FALSE), as.integer(0), PACKAGE="influenceR")
    return(x)
  }
  
  if ("package:Rmpi" %in% search()) {
    np <-  mpi.universe.size() - 1
    cl <- makeMPIcluster(np)
    
    x <- clusterApply(cl, 1:np, function(rank, el_i, n, m) {
        
        library(influenceR)
        .Call("snap_bridging_R", el_i, n, m, as.integer(TRUE), as.integer(rank), PACKAGE="influenceR")
        
    }, el_i, n, m) # ensure these values are exported.
    
    stopCluster(cl)
    mpi.exit()
    
    return(x[1])
  }
  else
    print("Error! Load Rmpi and supply a cluster object.")
    
}



ens <- function(g) {
  A <- get.adjacency(g)   # This will be sparse, which is great.
  S <- crossprod(A)       # S[i,j] = # of shared neighbors between i,j
  Q <- A * S              # Q[i,j] = # of shared neighbors if i and j are neighbors, 0 else
  qsum <- rowSums(Q)
  deg <- rowSums(A)
  ens <- deg - (qsum/deg)
}

constraint <- function(g) {
  
  process_sparse <- function(A, Ai, deg) {
    M <- as(A, 'TsparseMatrix')
    x <- .Call("process_sparse_R", M@i, M@j, M@x, Ai, deg, nnzero(M))
    M@x <- x
    M
  }
  
  A <- get.adjacency(g, sparse=T)
  n <- dim(A)[1]
  deg <- rowSums(A)

  constraint_i <- function(i) {
    
    #jqd <- (A * A[i,]) * deg
    #jqd[,i] <- 0
    #jqd@x <- 1/jqd@x  # In place reciprocal. Alternative: jqd[jqd==0] <- Inf; jqd <- 1/jqd
    
    jqd <- process_sparse(A, A[i,], deg)       
    Sj <- colSums(jqd)
  
    Cj <- Sj/deg[i] + (1/deg[i])
    Cj2 <- Cj * Cj
    sum(Cj2)
  }
  
  sapply(1:n, constraint_i)
}


main <- function(args) {
  
  graph <- dimacs.to.graph(args[[2]])
  func <- args[[1]]
  centrality <- c(eigen=eigencentrality, betweenness=betweenness, ens=ens, constraint=constraint, bridging=bridging)
  
  x <- centrality[[func]](graph)
  
  write.table(x, quote=F, row.names=F, col.names=F)
}
