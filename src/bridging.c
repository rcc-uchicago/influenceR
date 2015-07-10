#include <mpi.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

#include <graph_defs.h>
#include <graph_gen.h>


double bridging_vertex_precomp(graph_t *G, long v, double cls, double *closeness);
double *main_bridging(graph_t *G, int *edgelist, double *scores);
double closeness(graph_t *G, long ignore_edge0, long ignore_edge1);
long BFS_parallel_frontier_expansion_bridging(graph_t* G, long src, long diameter, double *distance, long ignore_edge0, long ignore_edge1 );


double *bridging(graph_t *G, int *edgelist, double *scores)
{  
  
 	int n = G->n; /* number of nodes */
	int m = G->m; /* number of edges */

  
  long u, v, j, k;
  
	/* 1) compute closeness by edge in file */
	
  double *closeness_by_edge = (double *) R_alloc(m, sizeof(double));
  
  for (int i = 0; i < m/2; i++) {
  
    u = edgelist[i*2] - 1;
    v = edgelist[i*2+1] - 1;
        
    /* Find edge numbers */
    for (j=G->numEdges[u]; v != G->endV[j] && j<G->numEdges[u+1]; j++);
    for (k=G->numEdges[v]; u != G->endV[k] && k<G->numEdges[v+1]; k++);
    assert(j != G->numEdges[u+1]);
    assert(k != G->numEdges[v+1]);
    
    /* Calculate closeness */
    double c = closeness(G, j, k);
    closeness_by_edge[j] = c;
    closeness_by_edge[k] = c;
  }

  /* 2) Compute closeness by vertex */
  
	double cls = closeness(G, -1, -1); // normal closeness (use all edges)
	for (v = 0; v < n; v++) 
		scores[v] = bridging_vertex_precomp(G, v, cls, closeness_by_edge);
  
  return scores;
}

double *bridging_MPI(graph_t *G, int *edgelist, double *scores)
{  
  
  // Get the number of processes
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //#ifdef VERBOSE
  fprintf(stderr, "hello from main_brdiging, process %d\n", rank);
  //#endif
  
 	int n = G->n; /* number of nodes */
	int m = G->m; /* number of edges */

  
	/* 1) compute closeness by edge in file */
	
  int bufsize = ceil(((double)m) / size), delta = bufsize/2;
  int start = rank * delta, end = start + delta;
  end = end > m/2 ? m/2 : end;
  
#ifdef VERBOSE
  fprintf(stderr, "%d range: %d-%d\n", rank, start, end); 
#endif
  
  double *buf = (double *) R_alloc(bufsize, sizeof(double));
  int *edgeidx = (int *) R_alloc(bufsize, sizeof(int));
  
  assert(buf);
  assert(edgeidx);
  
  int i=0, u, v;
  long j, k;
  
  
  for (int ii = start; ii < end; ii++) {
    u = edgelist[ii*2] - 1;
    v = edgelist[ii*2+1] - 1;
    
    /* Find edge numbers */
    for (j=G->numEdges[u]; v != G->endV[j] && j<G->numEdges[u+1]; j++);
    for (k=G->numEdges[v]; u != G->endV[k] && k<G->numEdges[v+1]; k++);
    assert(j != G->numEdges[u+1]);
    assert(k != G->numEdges[v+1]);
    
    /* Calculate closeness */
    buf[i] = closeness(G, j, k);
    edgeidx[i] = j;
    buf[i+1] = buf[i];
    edgeidx[i+1] = k;
    i+=2;

    //fprintf(stderr, "%d: CBE %d %d %g\n", rank, j, k, buf[i]);
  }

#ifdef VERBOSE
  fprintf(stderr, "Rank %d done reading edges\n", rank);
#endif
  

  double *closeness_buf = NULL;
  int *edge_indices = NULL;
  if (rank == 0) {
    closeness_buf = (double *) R_alloc(bufsize*size, sizeof(double));
    edge_indices = (int *) R_alloc(bufsize*size, sizeof(int));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Gather(buf, bufsize, MPI_DOUBLE, closeness_buf, bufsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(edgeidx, bufsize, MPI_INT, edge_indices, bufsize, MPI_INT, 0, MPI_COMM_WORLD);
  
  double *closeness_by_edge = (double *) R_alloc(m, sizeof(double));
  /* Fill REAL closeness_by_edge matrix */
    
  if (rank == 0) {
    for (int i = 0; i < m; i++) {
  	  closeness_by_edge[edge_indices[i]] = closeness_buf[i];
#ifdef VERBOSE
      printf("CBE %d %g\n", edge_indices[i], closeness_buf[i]); 
#endif
    }
    
    //free(closeness_buf);
    //free(edge_indices);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(closeness_by_edge, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //free(buf);
  //free(edgeidx);
	

	/* 2) compute bridging score by NODE. Parallization here may be more trouble than it's worth. But we already have the resources. */
  delta = ceil(((double)n) / size);
  start = rank * delta, end = start + delta;
  end = end > n ? n : end;
  
	double cls = closeness(G, -1, -1); // normal closeness (use all edges)

  buf = (double *) R_alloc(delta, sizeof(double));

	for (int v = start; v < end; v++) 
		buf[v-start] = bridging_vertex_precomp(G, v, cls, closeness_by_edge);
	
  MPI_Gather(buf, delta, MPI_DOUBLE, scores, delta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  return scores;
}

double bridging_vertex_precomp(graph_t *G, long v, double cls, double *closeness) {

  int n = G->n;

  int degree = 0;
  double sum = 0;

  for (long j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
    double cls_ = closeness[j];
  	sum += cls - cls_;
    degree++;
  }

  if (degree == 0)
  return 0;

  return sum/((double) degree);
}


// Two edges correspond to the same edge.
double closeness(graph_t *G, long ignore_edge0, long ignore_edge1)
{
	int n = G->n;
	
	double *distance = (double *) R_alloc(sizeof(double), n);
  assert(distance);
  if (distance == NULL) {
    fprintf(stderr, "RAN OUT OF MEM\n");
  }
  
	double sum = 0;
	
	for (int i = 0; i < n; i++) {
		/* memset */
		for (int j = 0; j < n; j++)
			distance[j] = INFINITY;
		
		BFS_parallel_frontier_expansion_bridging(G, i, 75, distance, ignore_edge0, ignore_edge1);
		
		for (int j = 0; j < i; j++) { /* sum upper triangular part */
		    sum += (1/distance[j]);
		}
	}

	//free(distance);
	return sum / (n*n - n);
}


/* 
 * OpenMP is disabled because we're parallelizing over graph edges and computing multiple BFS at once.
 * To enable, set the following macro directive:
 * DEFINE _OPENMP_BRIDGING _OPENMP
 * 
 */
long BFS_parallel_frontier_expansion_bridging(graph_t* G, long src, long diameter, double *distance, long ignore_edge0, long ignore_edge1 ) {

    attr_id_t* S;
    long *start;
    char* visited;
    long *pSCount;
#ifdef DIAGNOSTIC
    double elapsed_time;
#endif
#ifdef _OPENMP_BRIDGING
    omp_lock_t* vLock;
#endif

    long phase_num, numPhases;
    long count;


#ifdef _OPENMP_BRIDGING 

OMP("omp parallel")
    {
#endif

        attr_id_t *pS, *pSt;
        long pCount, pS_size;
        long v, w;
        int tid, nthreads;
        long start_iter, end_iter;    
        long j, k, vert, n;
#ifdef _OPENMP_BRIDGING
        int myLock;
#endif

#ifdef _OPENMP_BRIDGING    
        long i;
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif


#ifdef DIAGNOSTIC    
        if (tid == 0)
            elapsed_time = get_seconds();
#endif

        if (tid == 0)  
            numPhases = diameter + 1;
        n = G->n;

        pS_size = n/nthreads + 1;
        pS = (attr_id_t *) malloc(pS_size*sizeof(attr_id_t));
        assert(pS != NULL);

        if (tid == 0) {  
            S = (attr_id_t *) malloc(n*sizeof(attr_id_t));
            visited = (char *) calloc(n, sizeof(char));
            start = (long *) calloc((numPhases+2), sizeof(long));
            pSCount = (long *) malloc((nthreads+1)*sizeof(long));
#ifdef _OPENMP_BRIDGING
            vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
#endif
        }

#ifdef _OPENMP_BRIDGING    
OMP("omp barrier")
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_init_lock(&vLock[i]);
        }
#endif

#ifdef _OPENMP_BRIDGING
OMP("omp barrier")
#endif

        if (tid == 0) {
            S[0] = src;
            visited[src] = (char) 1;
            count = 1;
            phase_num = 0;
            start[0] = 0;
            start[1] = 1;
			distance[src] = 0;
        }


#ifdef _OPENMP_BRIDGING
OMP("omp barrier")
#endif

        while (start[phase_num+1] - start[phase_num] > 0) {

            pCount = 0;

            start_iter = start[phase_num];
            end_iter = start[phase_num+1];
#ifdef _OPENMP_BRIDGING
OMP("omp for")
#endif
            for (vert=start_iter; vert<end_iter; vert++) {

                v = S[vert];

                for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                   
                    if(j == ignore_edge0 || j == ignore_edge1) {
                        continue;
                    }

                    w = G->endV[j]; 
                    if (v == w)
                        continue;
#ifdef _OPENMP_BRIDGING
                    myLock = omp_test_lock(&vLock[w]);
                    if (myLock) {
#endif
                        if (visited[w] != (char) 1) { 
							distance[w] = distance[v] + 1;
                            visited[w] = (char) 1;
                            if (pCount == pS_size) {
                                /* Resize pS */
                                pSt = (attr_id_t *)
                                    malloc(2*pS_size*sizeof(attr_id_t));
                                memcpy(pSt, pS, pS_size*sizeof(attr_id_t));
                                free(pS);
                                pS = pSt;
                                pS_size = 2*pS_size;
                            }
                            pS[pCount++] = w;
                        }
#ifdef _OPENMP_BRIDGING
                        omp_unset_lock(&vLock[w]);
                    }
#endif
                }
            }


#ifdef _OPENMP_BRIDGING
OMP("omp barrier")
#endif            
            pSCount[tid+1] = pCount;

#ifdef _OPENMP_BRIDGING
OMP("omp barrier")
#endif            

            if (tid == 0) {
                pSCount[0] = start[phase_num+1];
                for(k=1; k<=nthreads; k++) {
                    pSCount[k] = pSCount[k-1] + pSCount[k];
                }
                start[phase_num+2] = pSCount[nthreads];
                count = pSCount[nthreads];
                phase_num++;
            }

#ifdef _OPENMP_BRIDGING
OMP("omp barrier")
#endif
            for (k = pSCount[tid]; k < pSCount[tid+1]; k++) {
                S[k] = pS[k-pSCount[tid]];
            } 


#ifdef _OPENMP_BRIDGING
OMP("omp barrier")
#endif
        } /* End of search */

#ifdef DIAGNOSTIC
        if (tid == 0) {
            fprintf(stderr, "Search from vertex %ld," 
                    " No. of vertices visited: %ld\n", src, count);
        }
#endif

        free(pS);
#ifdef _OPENMP_BRIDGING    
OMP("omp barrier")
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_destroy_lock(&vLock[i]);
        }
OMP("omp barrier")
#endif

        if (tid == 0) {
            free(S);
            free(start);
            free(visited);
            free(pSCount);
#ifdef _OPENMP_BRIDGING
            free(vLock);
#endif

        }

#ifdef _OPENMP_BRIDGING    
    }
#endif

#ifdef DIAGNOSTIC    
    elapsed_time = get_seconds() - elapsed_time;
    fprintf(stderr, "Time taken for BFS: %lf seconds\n", elpased_time);
#endif
    return count;
}



