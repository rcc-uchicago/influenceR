# igraph package has arpack functions
# evcent - eigenvector centrality
# betweenness - vertex betweenness

dimacs.to.graph <- function(fname) {
    x <- read.table(fname, skip=1)
    el <- as.matrix(x[2:3])
    
    graph.edgelist(el, directed=F)
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


bridging <- function(g, MPI=F, cluster=NULL) {
  
  el <- get.edgelist(g)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el_i))
  m <- as.integer(length(el_i)/2)
  
  if(!MPI) {
    x <- .Call("snap_bridging_R", el_i, n, m, as.integer(FALSE), as.integer(0), PACKAGE="influenceR")
    return(x)
  }
  
  if ("package:Rmpi" %in% search()) {
    x <- clusterCall(cl, function(...) {
        
        library(influenceR)
        .Call("snap_bridging_R", ..., PACKAGE="influenceR")
        
    }, el_i, n, m, as.integer(TRUE), as.integer(0)) # ensure these values are exported.
    
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
  centrality <- c(eigen=eigencentrality, betweenness=betweenness, ens=ens, constraint=constraint)
  
  x <- centrality[[func]](graph)
  
  write.table(x, quote=F, row.names=F, col.names=F)
}
