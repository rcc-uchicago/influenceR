#include <graph_defs.h>
#include <graph_gen.h>
#include <graph_metrics.h>

#include "bridging.h"
#include "keyplayer.h"

#include <R.h>
#include <Rinternals.h>

int read_graph_from_edgelist(graph_t* G, int *EL, long n, long m) {

    long i;
    long count, offset;
    int int_wt=1, *int_weight; // maybe some graphs have weights.
    long u, v;
    attr_id_t *src;
    attr_id_t *dest;
    attr_id_t *degree;

    // R will free this memory on return from .Call
    // see: http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Memory-allocation
    src = (attr_id_t *) R_alloc (m, sizeof(attr_id_t));
    dest = (attr_id_t *) R_alloc(m, sizeof(attr_id_t));
    degree = (attr_id_t *) R_alloc(n, sizeof(attr_id_t));
    for (int i = 0; i < n; i++) degree[i] = 0;

    assert(src != NULL);
    assert(dest != NULL);
    assert(degree != NULL);

    int_weight = (int *) malloc(m * sizeof(int));
    assert(int_weight != NULL);

    count = 0;

    for (int i = 0; i < m; i++) {
      u = EL[2*i];
      v = EL[2*i+1];
      
      if ((u <= 0) || (u > n) || (v <= 0) || (v > n)) {
          fprintf(stderr, "Error reading edge # %ld (%ld, %ld) in the input file."
                  " Please check the input graph file.\n", count+1, u, v);
          return 1;
      }
      src[count] = u-1;
      dest[count] = v-1;
      degree[u-1]++;
      degree[v-1]++;
      int_weight[count] = int_wt;
      
      count++;
    }

    if (count != m) {
        fprintf(stderr, "Error! Number of edges specified in problem line (%ld)" 
                " does not match the total number of edges (%ld) in file."
                " Please check"
                " the graph input file.\n", m, count);
        return 1;
    }

    /*
       for (i=0; i<m; i++) {
       fprintf(stderr, "[%d %d] ", src[i], dest[i]);
       }
     */

    G->endV = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
    assert(G->endV != NULL);

    G->edge_id = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
    assert(G->edge_id != NULL);

    G->numEdges = (attr_id_t *) malloc((n+1)*sizeof(attr_id_t));
    assert(G->numEdges != NULL);

    G->undirected = 1;
    G->weight_type = 1;
    G->zero_indexed = 0;

    G->n = n;
    G->m = 2*m;

    G->int_weight_e = (int *) malloc(G->m * sizeof(int));       
    assert(G->int_weight_e != NULL);

    /* ToDo: parallelize this step */
    G->numEdges[0] = 0; 
    for (i=1; i<=G->n; i++) {
        G->numEdges[i] = G->numEdges[i-1] + degree[i-1];
    }

    for (i=0; i<count; i++) {
        u = src[i];
        v = dest[i];
        
        offset = degree[u]--;
        G->endV[G->numEdges[u]+offset-1] = v;
        G->int_weight_e[G->numEdges[u]+offset-1] = int_weight[i];
        G->edge_id[G->numEdges[u]+offset-1] = i;

        offset = degree[v]--;
        G->endV[G->numEdges[v]+offset-1] = u;
        G->int_weight_e[G->numEdges[v]+offset-1] = int_weight[i];
        G->edge_id[G->numEdges[v]+offset-1] = i;
    } 

    /*          
    for (i=0; i<G->n; i++) {
        for (j=G->numEdges[i]; j<G->numEdges[i+1]; j++) {
            fprintf(stderr, "<%ld %ld %d> ", i+1, G->endV[j]+1, 
            G->int_weight_e[j]);
        }
    }
    */


    // No need to free code allocated with R_alloc
    //free(buf);
    //free(degree);
    //free(src);
    //free(dest);
    
    return 0;
}

int snap_betweenness(int *E, long n, long m, double *BC) {
  graph_t G;
  int r = read_graph_from_edgelist(&G, E, n, m);
  if (r) {
    printf("Error code!");
    return 1;
  }
  vertex_betweenness_centrality(&G, BC, n);
  return 0;
}

SEXP snap_betweenness_R(SEXP sE, SEXP sn, SEXP sm)
{
  
  int n = INTEGER(sn)[0],
    m = INTEGER(sm)[0];
  
  SEXP sBC = PROTECT(allocVector(REALSXP, n));
  
  int *E = INTEGER(sE);
  
  double *BC = REAL(sBC);
  
  snap_betweenness(E, n, m, BC);
  
  UNPROTECT(1);
 
  return sBC;
}

SEXP snap_keyplayer_R(SEXP sE, SEXP sn, SEXP sm, SEXP sk, SEXP sprob, SEXP stol, SEXP sMaxsec, SEXP sRoundsec) {
  int n = INTEGER(sn)[0],
    m = INTEGER(sm)[0],
    k = INTEGER(sk)[0];
    
  double prob = REAL(sprob)[0],
    tol = REAL(stol)[0];
    
  long maxsec = INTEGER(sMaxsec)[0],
    roundsec = INTEGER(sRoundsec)[0];
  
  int *E = INTEGER(sE);
  graph_t G;
  int r = read_graph_from_edgelist(&G, E, n, m);
  
  
  SEXP sKP = PROTECT(allocVector(INTSXP, n));
  int *KP = INTEGER(sKP);

#ifdef OPENMP
  keyplayer_driver_omp(&G, n, k, prob, tol, maxsec, roundsec, KP);
#else
  keyplayer_driver(&G, n, k, prob, tol, maxsec, KP);
#endif
  
  
  UNPROTECT(1);
 
  return sKP;
}

/* Rank is rank if this is an MPI call, 0 else */
SEXP snap_bridging_R(SEXP sE, SEXP sn, SEXP sm, SEXP sMPI, SEXP srank) {
  int n = INTEGER(sn)[0],
    m = INTEGER(sm)[0],
    rank = INTEGER(srank)[0],
    mpi = INTEGER(sMPI)[0]; // use mpi?
  
  int *E = INTEGER(sE);
  graph_t G;
  int r = read_graph_from_edgelist(&G, E, n, m);
  
  SEXP sBC = PROTECT(allocVector(REALSXP, rank==0 ? n : 0));
  if (rank == 0) {
    //sBC = PROTECT(allocVector(REALSXP, n));
#ifdef VERBOSE
    printf("Rank %d: allocated memory for %d doubles\n", rank, n);
#endif
    if (REAL(sBC) == NULL) {
        printf("Rank %d: error!\n", rank);
        exit(1);
    }
  }
  else {
    printf("Rank %d: Did not allocate memory\n", rank);
  }
  double *BC = REAL(sBC);

#ifdef USE_MPI
  if (mpi != 0)
    bridging_MPI(&G, E, BC);
  else
    bridging(&G, E, BC);
#else
  bridging(&G, E, BC);
#endif
  
  UNPROTECT(1);
 
  return sBC;
}
