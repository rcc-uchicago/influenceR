# graph_metrics.R: R code for vertex importance metrics.
# AUTHOR: Simon Jacobs <sdjacobs@uchicago.edu>
# LICENSE: GPLv2


#' Convert a CSV file to an igraph graph object.
#'
#' The first column should be sources, the second should be targets.
#'
#' @param fname A filename
#' @return An igraph graph object built from the filename.
#' 
#' @export
csv.to.igraph <- function(fname) {
    x <- read.csv(fname) # this may be dangerous because of users' settings.
                         # See: http://r-pkgs.had.co.nz/r.html
    el <- as.matrix(x[c(1,2)])
    if(!is.character(el))
      el <- apply(el, 2, as.character)
    
    igraph::graph.edgelist(el, directed=F)
}

#' Vertex betweenness centrality measure.
#'
#' The betweenness centrality score of a node u is the sum over all pairs s,t of the
#' proportion of shortest paths between s and t that pass through u. This 
#' function allows the use of either the SNAP betweenness implementation (default), or 
#' the igraph betweenness function. The SNAP version makes use of OpenMP for 
#' parallelization, and may be faster in some circumstances.
#'
#' @seealso \url{http://snap-graph.sourceforge.net/}
#'
#' @param g The igraph object to analyze
#' @param snap True to use the SNAP betweenness code, False to use igraph::betweenness
#' @return A numeric vector with the betweenness centrality score for each vertex
#'
#' @export
betweenness <- function(g, snap=T) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
# 1/2 the values of our betweenness code, which is because this is UNDIRECTED for real
  if (!snap)
    return(igraph::betweenness(g))
  
  el <- igraph::get.edgelist(g, names=F)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el))
  m <- as.integer(length(el)/2) # TODO: for directed too?
  
  vals <- .Call("snap_betweenness_R", el_i, n, m, PACKAGE="influenceR")
  vals[vals<2^-128] <- 0
  vals[is.nan(vals)] <- 0
  names(vals) <- igraph::V(g)$name
  vals
}

#' Compute a KPP-Pos set for a given graph.
#'
#' The "Key Player" family of node importance algorithms (Borgatti 2006) involves the selection
#' of a metric of node importance and a combinatorial optimization strategy to
#' choose the set S of vertices of size k that maximize that metric. This
#' function implements KPP-Pos, an algorithm to to identify |S| actors that optimize information diffusion
#' through the network. We sum over all vertices not in S the reciprocal
#' of the shortest distance to a vertex in S. For combinatorial optimization, we use
#' stochastic gradient descent, where in each optimization round, we select a node u in S
#' and v not in S at random, switch them, and accept the switch if evaluation of the
#' metric improves. This implementation uses OpenMP (if available on the host system) so that
#' multiple workers can explore the solution space in parallel, synchronizing to pick the best
#' answer after a given computation budget has elapsed.
#'
#' @seealso \url{http://www.bebr.ufl.edu/sites/default/files/Borgatti\%20-\%202006\%20-\%20Identifying\%20sets\%20of\%20key\%20players\%20in\%20a\%20social\%20networ.pdf}
#'
#' @param g The igraph object to analyze.
#' @param k The size of the KP-set
#' @param prob probability of accepting a state with a lower value
#' @param tol tolerance within which to stop the optimization and accept the current value
#' @param maxsec The total computation budget for the optimization, in seconds
#' @param roundsec Number of seconds in between synchronizing workers' answer
#' @return a vector with the vertex number of each vertex in the selected set S.
#'
#' @export
keyplayer <- function(g, k, prob = 0.0, tol = 0.0001, maxsec = 600, roundsec = 30) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  el <- igraph::get.edgelist(g, names=F)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el))
  m <- as.integer(length(el)/2)

  s <- .Call("snap_keyplayer_R", el_i, n, m, as.integer(k), prob, tol, as.integer(maxsec), as.integer(roundsec), PACKAGE="influenceR")
  
  igraph::V(g)[which(s>0)]
}

#' Valente's bridging vertex measure.
#'
#' A node's bridging score is the average decrease in cohesiveness if each of
#' its edges were removed from the graph.
#' 
#' @seealso \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2889704/}
#'
#' @param g The igraph object to analyze.
#' @return A numeric vector with the bridging score for each vertex
#'
#' @export
bridging <- function(g) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  el <- igraph::get.edgelist(g, names = F)
  el_i <- as.integer(t(el))
  n <- as.integer(max(el_i))
  m <- as.integer(length(el_i) / 2)
  
  x <- .Call("snap_bridging_R", el_i, n, m, as.integer(FALSE), as.integer(0), PACKAGE = "influenceR")
  names(x) <- igraph::V(g)$name
  x
}

#' Burt's effective network size vertex measure.
#'
#' @param g The igraph object to analyze.
#' @return A numeric vector with the effective network size for each vertex
#'
#' @export
ens <- function(g) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  A <- igraph::get.adjacency(g)   # This will be sparse, which is great.
  S <- Matrix::crossprod(A)       # S[i,j] = # of shared neighbors between i,j
  Q <- A * S              # Q[i,j] = # of shared neighbors if i and j are neighbors, 0 else
  qsum <- Matrix::rowSums(Q)
  deg <- Matrix::rowSums(A)
  ens <- deg - (qsum / deg)
  ens[is.nan(ens)] <- 0 # If a vertex has no neighbors, make its ENS 0
  names(ens) <- igraph::V(g)$name
  ens
}

#' Burt's graph constraint vertex measure.
#'
#' This is an alternative to the implementation of Burt's consraint in the
#' igraph package.
#'
#' @param g The igraph object to analyze.
#' @param v vertices over which to compute constraint (default to all)
#' @return A numeric vector with the constraint score for each vertex in v
#'
#' @export
constraint <- function(g, v=igraph::V(g)) {
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  
  process_sparse <- function(A, Ai, deg) {
    M <- as(A, 'TsparseMatrix')
    x <- .Call("process_sparse_R", M@i, M@j, M@x, Ai, deg, Matrix::nnzero(M), PACKAGE = "influenceR")
    M@x <- x
    M
  }
  
  A <- igraph::get.adjacency(g, sparse=T)
  n <- dim(A)[1]
  deg <- Matrix::rowSums(A)

  constraint_i <- function(i) {
    # process sparse does this: jq <- drop0(t(A*A[,i]) * A[,i]); jqd <- drop0(jq * deg)
    jqd <- process_sparse(A, A[i, ], deg)
    
    jqd <- Matrix::drop0(jqd)
    jqd@x <- (1 / jqd@x) * (1 / deg[i])
         
    Sj <- Matrix::colSums(jqd)
  
    idx <- as.numeric(igraph::neighbors(g, i))
    Sj[idx] <- Sj[idx] + (1 / deg[i])
    
    Sj2 <- Sj * Sj
    sum(Sj2)
  }
  
  vals <- sapply(v, constraint_i)
  names(vals) <- v$name
  vals
}
