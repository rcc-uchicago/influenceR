#include "graph_defs.h"
#include "graph_kernels.h"

long graph_diameter(graph_t* G) {

#if 0
    long i, j, n;
    long tid, nthreads;
    long* S;          /* The BFS queue */ 
    long* d;
    long start, end;
    long *local_dia, dia;
    long v, w;
    double running_time;
    /* The outer loop is parallelized in this case. Each thread does a BFS 
    from a vertex, and the closeness centrality value is computed */   
#ifdef _OPENMP
OMP("omp parallel "
    "private(tid, nthreads, S, d, start, end, v, w) "
    "private(i, j, n, running_time)")
{
#endif

#ifdef _OPENMP	
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads(); 
#else
	tid = 0;
	nthreads= 1;
#endif

	n = G->n;

    S   = (long *) malloc(n*sizeof(long));
    assert(S != NULL);
  
    d = (long *) malloc(n*sizeof(long));
    assert(d != NULL);
   
    for (i=0; i<n; i++) {
        d[i] = -1;    
    }
 
    if (tid == 0) {
        local_dia = (long *) calloc(nthreads, sizeof(long));
        dia = 0;
#ifdef _OPENMP
		running_time = omp_get_wtime();
#endif
	}

#ifdef _OPENMP	
    #pragma omp barrier
#endif

#ifdef _OPENMP	
    #pragma omp for schedule(dynamic)
#endif
	for (i=0; i<n; i++) {
#if VERBOSE_OUTPUT
        if ((i % 1000) == 0)
            fprintf(stdout, "%ld ", i);
#endif
        start = 0;
        S[0] = i;
        end = 1;
        d[i] = 0;

        while (end - start > 0) {
            v = S[start];
            for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                w = G->endV[j];
                if (v != w) {
                    /* w found for the first time? */ 
                    if (d[w] < 0) {
                        S[end++] = w;
                        d[w] = d[v] + 1;
                    }
                }
            }
            start++;
        }
       
        if (d[S[end-1]] > local_dia[tid])
           local_dia[tid] = d[S[end-1]];

        for (j=0; j<end; j++) {
            w = S[j];
            d[w] = -1;
        }
    }

#ifdef _OPENMP	
    #pragma omp barrier
#endif
    if (tid == 0) {
        dia = local_dia[0];
        for (j=1; j<nthreads; j++) {
            if (local_dia[j] > dia)
                dia = local_dia[j];
        }
        *diaPtr = dia;
    }

#ifdef _OPENMP
    #pragma omp barrier
#endif
	
    if (tid == 0) { 
#ifdef _OPENMP
        running_time = omp_get_wtime() - running_time;
        fprintf(stderr, "done.\nTime taken: %lf seconds.\n", running_time);
#endif
		free(local_dia);
    }
    
    free(S);
    free(d);
#ifdef _OPENMP
	#pragma omp barrier
#endif

#ifdef _OPENMP	
}    
#endif
#endif
    return 0;
}
