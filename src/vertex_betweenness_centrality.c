#include "graph_defs.h"
#include "graph_metrics.h"
#include "prefix_sums.h"

/* We replace sprng from the SNAP code with R RNG: */
#include <R.h>
#define SPRNG_DEFAULT 0 
#define init_sprng(a, b, c, d, e) NULL; omp_set_lock(&rnglock); GetRNGstate(); omp_unset_lock(&rnglock)
#define sprng(a) unif_rand()
#define free_sprng(a) omp_set_lock(&rnglock); PutRNGstate(); omp_unset_lock(&rnglock)
#define fprintf(a, ...) REprintf(__VA_ARGS__)
#define exit(x) error("SNAP code exited with error: %d\n", x)

void vertex_betweenness_centrality_parBFS(graph_t* G, double* BC, long numSrcs) {


    omp_lock_t rnglock;
    omp_init_lock(&rnglock);

    attr_id_t *S;      /* stack of vertices in the order of non-decreasing 
                          distance from s. Also used to implicitly 
                          represent the BFS queue */
    plist_t* P;        /* predecessors of a vertex v on shortest paths from s */
    double* sig;       /* No. of shortest paths */
    attr_id_t* d;      /* Length of the shortest path between every pair */
    double* del;       /* dependency of vertices */
    attr_id_t *in_degree, *numEdges, *pSums;
    attr_id_t* pListMem;    
#if RANDSRCS
    attr_id_t* Srcs; 
#endif
    attr_id_t *start, *end;
    long MAX_NUM_PHASES;
    attr_id_t *psCount;

#ifdef _OPENMP    
    omp_lock_t* vLock;
    long chunkSize;
#endif
#ifdef DIAGNOSTIC
    double elapsed_time;
#endif
    int seed = 2387;

#ifdef _OPENMP    
OMP("omp parallel firstprivate(G)")
    {
#endif

        attr_id_t *myS, *myS_t;
        attr_id_t myS_size;
        long i, j, k, p, count, myCount;
        long v, w, vert;
        long k0, k1;
        long numV, num_traversals, n, m, phase_num;
        long start_iter, end_iter;
        long tid, nthreads;
        int* stream;
#ifdef DIAGNOSTIC
        double elapsed_time_part;
#endif

#ifdef _OPENMP
        int myLock;
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time = get_seconds();
            elapsed_time_part = get_seconds();
        }
#endif

        /* numV: no. of vertices to run BFS from = numSrcs */
        numV = numSrcs;
        n = G->n;
        m = G->m;

        /* Permute vertices */
        if (tid == 0) {
#if RANDSRCS
            Srcs = (attr_id_t *) malloc(n*sizeof(attr_id_t));
#endif
#ifdef _OPENMP
            vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
#endif
        }

#ifdef _OPENMP   
OMP("omp barrier")
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_init_lock(&vLock[i]);
        }
#endif

        /* Initialize RNG stream */ 
        stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);

#if RANDSRCS
#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            Srcs[i] = i;
        }

#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            j = n * sprng(stream);
            if (i != j) {
#ifdef _OPENMP
                int l1 = omp_test_lock(&vLock[i]);
                if (l1) {
                    int l2 = omp_test_lock(&vLock[j]);
                    if (l2) {
#endif
                        k = Srcs[i];
                        Srcs[i] = Srcs[j];
                        Srcs[j] = k;
#ifdef _OPENMP  
                        omp_unset_lock(&vLock[j]);
                    }
                    omp_unset_lock(&vLock[i]);
                }
#endif        
            }
        } 
#endif

#ifdef _OPENMP    
OMP("omp barrier")
#endif

        if (tid == 0) {
            MAX_NUM_PHASES = 500;
        }

#ifdef _OPENMP
OMP("omp barrier")
#endif

        /* Initialize predecessor lists */

        /* The size of the predecessor list of each vertex is bounded by 
           its in-degree. So we first compute the in-degree of every
           vertex */ 

        if (tid == 0) {
            P   = (plist_t  *) calloc(n, sizeof(plist_t));
            in_degree = (attr_id_t *) calloc(n+1, sizeof(attr_id_t));
            numEdges = (attr_id_t *) malloc((n+1)*sizeof(attr_id_t));
            pSums = (attr_id_t *) malloc(nthreads*sizeof(attr_id_t));
        }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
        for (i=0; i<m; i++) {
            v = G->endV[i];
#ifdef _OPENMP
            omp_set_lock(&vLock[v]);
#endif
            in_degree[v]++;
#ifdef _OPENMP
            omp_unset_lock(&vLock[v]);
#endif
        }

        prefix_sums(in_degree, numEdges, pSums, n);

        if (tid == 0) {
            pListMem = (attr_id_t *) malloc(m*sizeof(attr_id_t));
        }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            P[i].list = pListMem + numEdges[i];
            P[i].degree = in_degree[i];
            P[i].count = 0;
        }

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() -elapsed_time_part;
            fprintf(stderr, "In-degree computation time: %lf seconds\n", 
                    elapsed_time_part);
            elapsed_time_part = get_seconds();
        }
#endif

        /* Allocate shared memory */ 
        if (tid == 0) {
            free(in_degree);
            free(numEdges);
            free(pSums);

            S   = (attr_id_t *) malloc(n*sizeof(attr_id_t));
            sig = (double *) malloc(n*sizeof(double));
            d   = (attr_id_t *) malloc(n*sizeof(attr_id_t));
            del = (double *) calloc(n, sizeof(double));

            start = (attr_id_t *) malloc(MAX_NUM_PHASES*sizeof(attr_id_t));
            end = (attr_id_t *) malloc(MAX_NUM_PHASES*sizeof(attr_id_t));
            psCount = (attr_id_t *) malloc((nthreads+1)*sizeof(attr_id_t));
        }

        /* local memory for each thread */  
        myS_size = (2*n)/nthreads;
        myS = (attr_id_t *) malloc(myS_size*sizeof(attr_id_t));
        num_traversals = 0;
        myCount = 0;

#ifdef _OPENMP    
OMP("omp barrier")
#endif

#ifdef _OPENMP    
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            d[i] = -1;
        }

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stderr, "BC initialization time: %lf seconds\n", 
                    elapsed_time_part);
            elapsed_time_part = get_seconds();
        }
#endif

        for (p=0; p<n; p++) {
#if RANDSRCS
            i = Srcs[p];
#else
            i = p;
#endif
            if (G->numEdges[i+1] - G->numEdges[i] == 0) {
                continue;
            } else {
                num_traversals++;
            }

            if (num_traversals == numV + 1) {
                break;
            }

            if (tid == 0) {
                sig[i] = 1;
                d[i] = 0;
                S[0] = i;
                start[0] = 0;
                end[0] = 1;
            }

            count = 1;
            phase_num = 0;

#ifdef _OPENMP       
OMP("omp barrier")
#endif

            while (end[phase_num] - start[phase_num] > 0) {

                myCount = 0;
                start_iter = start[phase_num];
                end_iter = end[phase_num];
#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for schedule(dynamic) nowait")
#endif
                for (vert = start_iter; vert < end_iter; vert++) {
                    v = S[vert];
                    for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {

                        w = G->endV[j];
                        if (v != w) {

#ifdef _OPENMP                            
                            myLock = omp_test_lock(&vLock[w]);
                            if (myLock) { 
#endif             
                                /* w found for the first time? */ 
                                if (d[w] == -1) {
                                    if (myS_size == myCount) {
                                        /* Resize myS */
                                        myS_t = (attr_id_t *)
                                            malloc(2*myS_size*sizeof(attr_id_t));
                                        memcpy(myS_t, myS, 
                                                myS_size*sizeof(attr_id_t));
                                        free(myS);
                                        myS = myS_t;
                                        myS_size = 2*myS_size;
                                    }
                                    myS[myCount++] = w;
                                    d[w] = d[v] + 1;
                                    sig[w] = sig[v];
                                    P[w].list[P[w].count++] = v;
                                } else if (d[w] == d[v] + 1) {
                                    sig[w] += sig[v];
                                    P[w].list[P[w].count++] = v;
                                }
#ifdef _OPENMP  

                                omp_unset_lock(&vLock[w]);
                            } else {
                                if ((d[w] == -1) || (d[w] == d[v]+ 1)) {
                                    omp_set_lock(&vLock[w]);
                                    sig[w] += sig[v];
                                    P[w].list[P[w].count++] = v;
                                    omp_unset_lock(&vLock[w]);
                                }
                            }
#endif

                        }
                    }
                }
                /* Merge all local stacks for next iteration */
                phase_num++; 
                if (tid == 0) {
                    if (phase_num >= MAX_NUM_PHASES) {
                        fprintf(stderr, "Error: Max num phases set to %ld\n",
                                MAX_NUM_PHASES);
                        fprintf(stderr, "Diameter of input network greater than"
                                " this value. Increase MAX_NUM_PHASES"
                                " in vertex_betweenness_centrality_parBFS()\n");
                        exit(-1);
                    }
                }

                psCount[tid+1] = myCount;

#ifdef _OPENMP
OMP("omp barrier")
#endif

                if (tid == 0) {
                    start[phase_num] = end[phase_num-1];
                    psCount[0] = start[phase_num];
                    for(k=1; k<=nthreads; k++) {
                        psCount[k] = psCount[k-1] + psCount[k];
                    }
                    end[phase_num] = psCount[nthreads];
                }



#ifdef _OPENMP
OMP("omp barrier")
#endif

                k0 = psCount[tid]; 
                k1 = psCount[tid+1];
                for (k = k0; k < k1; k++) {
                    S[k] = myS[k-k0];
                } 

                count = end[phase_num];
            }

            phase_num--;

            while (phase_num > 0) {
                start_iter = start[phase_num];
                end_iter = end[phase_num];
#ifdef _OPENMP        
OMP("omp for schedule(static) nowait")
#endif
                for (j=start_iter; j<end_iter; j++) {
                    w = S[j];
                    for (k = 0; k<P[w].count; k++) {
                        v = P[w].list[k];
#ifdef _OPENMP
                        omp_set_lock(&vLock[v]);
#endif
                        del[v] = del[v] + sig[v]*(1+del[w])/sig[w];
#ifdef _OPENMP
                        omp_unset_lock(&vLock[v]);
#endif
                    }
                    BC[w] += del[w];
                }

                phase_num--;

#ifdef _OPENMP
OMP("omp barrier")
#endif            
            }


#ifdef _OPENMP
            chunkSize = n/nthreads;
OMP("omp for schedule(static, chunkSize) nowait")
#endif
            for (j=0; j<count; j++) {
                w = S[j];
                d[w] = -1;
                del[w] = 0;
                P[w].count = 0;
            }


#ifdef _OPENMP
OMP("omp barrier")
#endif

        }

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stderr, "BC computation time: %lf seconds\n", 
                    elapsed_time_part);
        }
#endif


#ifdef _OPENMP
OMP("omp barrier")
#endif

#ifdef _OPENMP
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_destroy_lock(&vLock[i]);
        }
#endif

        free(myS);

        if (tid == 0) { 
            free(S);
            free(pListMem);
            free(P);
            free(sig);
            free(d);
            free(del);
#ifdef _OPENMP
            free(vLock);
#endif
            free(start);
            free(end);
            free(psCount);

#ifdef DIAGNOSTIC
            elapsed_time = get_seconds() - elapsed_time;
            fprintf(stderr, "Time taken: %lf\n seconds", elapsed_time);
#endif

#if RANDSRCS
            free(Srcs);
#endif
        }

        free_sprng(stream);
#ifdef _OPENMP
    }    
#endif

}

void vertex_betweenness_centrality_simple(graph_t* G, double* BC, long numSrcs) {

    omp_lock_t rnglock;
    omp_init_lock(&rnglock);

    attr_id_t *in_degree, *numEdges, *pSums;
#if RANDSRCS
    attr_id_t* Srcs; 
#endif
    long num_traversals = 0;
#ifdef _OPENMP    
    omp_lock_t* vLock;
    long chunkSize;
#endif
#ifdef DIAGNOSTIC
    double elapsed_time;
#endif
    int seed = 2387;

    /* The outer loop is parallelized in this case. Each thread does a BFS 
       and the vertex BC values are incremented atomically */   
#ifdef _OPENMP
OMP("omp parallel firstprivate(G)")
    {
#endif
        attr_id_t *S;      /* stack of vertices in the order of non-decreasing 
                              distance from s. Also used to implicitly 
                              represent the BFS queue */
        plist_t* P;          /* predecessors of a vertex v on shortest paths 
                                from s */
        attr_id_t* pListMem;    
        double* sig;       /* No. of shortest paths */
        attr_id_t* d;      /* Length of the shortest path between every pair */
        double* del;       /* dependency of vertices */
        attr_id_t *start, *end;
        long MAX_NUM_PHASES;

        long i, j, k, p, count;
        long v, w, vert;
        long numV, n, m, phase_num;
        long tid, nthreads;
        int* stream;
#ifdef DIAGNOSTIC
        double elapsed_time_part;
#endif

#ifdef _OPENMP
        int myLock;
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time = get_seconds();
            elapsed_time_part = get_seconds();
        }
#endif

        /* numV: no. of vertices to run BFS from = numSrcs */
        numV = numSrcs;
        n = G->n;
        m = G->m;

        /* Permute vertices */
        if (tid == 0) {
#if RANDSRCS
            Srcs = (attr_id_t *) malloc(n*sizeof(attr_id_t));
#endif
#ifdef _OPENMP
            vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
#endif
        }

#ifdef _OPENMP   
OMP("omp barrier")
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_init_lock(&vLock[i]);
        }
#endif

        /* Initialize RNG stream */ 
        stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);

#if RANDSRCS
#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            Srcs[i] = i;
        }

#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            j = n * sprng(stream);
            if (i != j) {
#ifdef _OPENMP
                int l1 = omp_test_lock(&vLock[i]);
                if (l1) {
                    int l2 = omp_test_lock(&vLock[j]);
                    if (l2) {
#endif
                        k = Srcs[i];
                        Srcs[i] = Srcs[j];
                        Srcs[j] = k;
#ifdef _OPENMP  
                        omp_unset_lock(&vLock[j]);
                    }
                    omp_unset_lock(&vLock[i]);
                }
#endif        
            }
        } 
#endif

#ifdef _OPENMP    
OMP("omp barrier")
#endif

        MAX_NUM_PHASES = 50;

        /* Initialize predecessor lists */

        /* The size of the predecessor list of each vertex is bounded by 
           its in-degree. So we first compute the in-degree of every
           vertex */ 

        if (tid == 0) {
            in_degree = (attr_id_t *) calloc(n+1, sizeof(attr_id_t));
            numEdges = (attr_id_t *) malloc((n+1)*sizeof(attr_id_t));
            pSums = (attr_id_t *) malloc(nthreads*sizeof(attr_id_t));
        }


#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
        for (i=0; i<m; i++) {
            v = G->endV[i];
#ifdef _OPENMP
            omp_set_lock(&vLock[v]);
#endif
            in_degree[v]++;
#ifdef _OPENMP
            omp_unset_lock(&vLock[v]);
#endif
        }

        prefix_sums(in_degree, numEdges, pSums, n);

        P  = (plist_t  *) calloc(n, sizeof(plist_t));
        pListMem = (attr_id_t *) malloc(m*sizeof(attr_id_t));

        for (i=0; i<n; i++) {
            P[i].list = pListMem + numEdges[i];
            P[i].degree = in_degree[i];
            P[i].count = 0;
        }

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() -elapsed_time_part;
            fprintf(stderr, "In-degree computation time: %lf seconds\n", 
                    elapsed_time_part);
            elapsed_time_part = get_seconds();
        }
#endif

#ifdef _OPENMP
OMP("omp barrier")
#endif

        /* Allocate shared memory */ 
        if (tid == 0) {
            free(in_degree);
            free(numEdges);
            free(pSums);
        }

        S   = (attr_id_t *) malloc(n*sizeof(attr_id_t));
        sig = (double *) malloc(n*sizeof(double));
        d   = (attr_id_t *) malloc(n*sizeof(attr_id_t));
        del = (double *) calloc(n, sizeof(double));

        start = (attr_id_t *) malloc(MAX_NUM_PHASES*sizeof(attr_id_t));
        end = (attr_id_t *) malloc(MAX_NUM_PHASES*sizeof(attr_id_t));

#ifdef _OPENMP   
OMP("omp barrier")
#endif

        for (i=0; i<n; i++) {
            d[i] = -1;
        }

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stderr, "BC initialization time: %lf seconds\n", 
                    elapsed_time_part);
            elapsed_time_part = get_seconds();
        }
#endif

#ifdef _OPENMP
OMP("omp for reduction(+:num_traversals)")
#endif
        for (p=0; p<numV; p++) {
#if RANDSRCS
            i = Srcs[p];
#else
            i = p;
#endif
            if (G->numEdges[i+1] - G->numEdges[i] == 0) {
                continue;
            } else {
                num_traversals++;
            }

            sig[i] = 1;
            d[i] = 0;
            S[0] = i;
            start[0] = 0;
            end[0] = 1;

            count = 1;
            phase_num = 0;

            while (end[phase_num] - start[phase_num] > 0) {

                for (vert = start[phase_num]; vert < end[phase_num]; vert++) {
                    v = S[vert];
                    for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                        w = G->endV[j];
                        if (v != w) {
                            /* w found for the first time? */ 
                            if (d[w] == -1) {
                                S[count++] = w;
                                d[w] = d[v] + 1;
                                sig[w] = sig[v];
                                P[w].list[P[w].count++] = v;
                            } else if (d[w] == d[v] + 1) {
                                sig[w] += sig[v];
                                P[w].list[P[w].count++] = v;
                            }
                        }
                    }
                }

                phase_num++; 

                start[phase_num] = end[phase_num-1];
                end[phase_num] = count;
            }

            phase_num--;

            while (phase_num > 0) {
                for (j=start[phase_num]; j<end[phase_num]; j++) {
                    w = S[j];
                    for (k = 0; k<P[w].count; k++) {
                        v = P[w].list[k];
                        del[v] = del[v] + sig[v]*(1+del[w])/sig[w];
                    }
#ifdef _OPENMP
                    omp_set_lock(&vLock[w]);
                    BC[w] += del[w];
                    omp_unset_lock(&vLock[w]);
#else
                    BC[w] += del[w];
#endif
                }

                phase_num--;
            }

            for (j=0; j<count; j++) {
                w = S[j];
                d[w] = -1;
                del[w] = 0;
                P[w].count = 0;
            }

        }

#ifdef DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stderr, "BC computation time: %lf seconds\n", 
                    elapsed_time_part);
        }
#endif


#ifdef _OPENMP
OMP("omp barrier")
#endif

#ifdef _OPENMP
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_destroy_lock(&vLock[i]);
        }
#endif

        free(S);
        free(pListMem);
        free(P);
        free(sig);
        free(d);
        free(del);
        free(start);
        free(end);

        if (tid == 0) {

#ifdef _OPENMP
            free(vLock);
#endif

#if RANDSRCS
            free(Srcs);
#endif

#ifdef DIAGNOSTIC
            elapsed_time = get_seconds() - elapsed_time;
            fprintf(stderr, "Total time taken: %lf seconds\n", elapsed_time);
#endif

        }

        free_sprng(stream);

#ifdef _OPENMP
OMP("omp barrier")
    }
#endif

}    

void vertex_betweenness_centrality(graph_t* G, double* BC, long numSrcs) {

#ifdef _OPENMP
    /* Exact BC */
    if (G->n == numSrcs) {
        if (G->n < 5000) {
            vertex_betweenness_centrality_simple(G, BC, numSrcs);    
        } else {
            vertex_betweenness_centrality_parBFS(G, BC, numSrcs);
        }
    } else if (numSrcs < 50) {
        vertex_betweenness_centrality_simple(G, BC, numSrcs);
    } else {
        vertex_betweenness_centrality_parBFS(G, BC, numSrcs);
    }
#else
    vertex_betweenness_centrality_simple(G, BC, numSrcs);
#endif

}
