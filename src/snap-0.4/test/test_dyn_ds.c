#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"
#include "sprng.h"

int main(int argc, char** argv) {

    dyn_graph_t* G;
    double *params;
    char* filename;
    attr_id_t *src, *dest, *degree;
    attr_id_t *del_src, *del_dest;
    attr_id_t num_dels, tp;
    attr_id_t *taken;
    attr_id_t ts1, ts2;

    if (argc != 2) {
        fprintf(stderr, "Usage: ./test_dyn_ds <config_file_name>\n");
        exit(-1);
    }

    filename = (char *) malloc(1000*sizeof(char)); 
    strcpy(filename, argv[1]);

    G = (dyn_graph_t *) calloc(1, sizeof(dyn_graph_t));
    params = (double *) malloc(10*sizeof(double));
    
    read_dyn_test_config_file(G, filename, params);
    
    src = (attr_id_t *) malloc (G->m * sizeof(attr_id_t));
    dest = (attr_id_t *) malloc(G->m * sizeof(attr_id_t));
    degree = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    num_dels = 10000;
    del_src = (attr_id_t *) malloc(num_dels * sizeof(attr_id_t));
    del_dest = (attr_id_t *) malloc(num_dels * sizeof(attr_id_t));
    taken = (attr_id_t *) calloc(G->m , sizeof(attr_id_t));

    assert(src != NULL);
    assert(dest != NULL);
    assert(degree != NULL);

    /* Generate edges */
    par_gen_RMAT_edges(G, params, src, dest, degree);

    /* Insertions only */

    /* Exp 1: Oracle, we know bucket sizes */
    // dyn_ds_init_nomalloc(G, src, dest, degree);
    /* Exp 2: We resize buckets inside code, locking */
    dyn_ds_init_malloc(G, src, dest);

    /* Exp 3: Vertex partitioning, no lock, no malloc */ 
    dyn_vpart_init_nomalloc(G, src, dest, degree);
    
    /* Exp 4: Edge partitioning, no lock, malloc */
    dyn_epart_init_nomalloc(G, src, dest, degree);

    /* Exp 5: Binary heaps, lock, no malloc */

    /* Exp 6: Sorting + after sorting assignment */
    // dyn_ds_init_base(G, src, dest, degree);

    /* Exp 7: Binary heap for high degree + arrays for low degree */
    // dyn_hybrid_init_nomalloc(G, src, dest, degree);
    
    // dyn_heap_init_nomalloc(G, src, dest, degree);
    /* Deletions */

    /* Create deletion tuples */
    /*
    tp = 0;   

    while (tp < num_dels) {
        int k = (int) (drand48() * G->m);
        if (taken[k] == 0) {
            del_src[tp] = src[k];
            del_dest[tp] = dest[k];
            taken[k] = 1;
            tp++;
        }
        
    }
    */
    /* Exp 8: Dyn array deletions */
    // dyn_gr_init(G, src, dest, degree);
    // dyn_gr_del(G, del_src, del_dest, degree, num_dels);
    
    /* Exp 9: Heap deletions */
    // dyn_heap_del(G, src, dest, num_del);

    /* Exp 10: Hybrid DS deletions */
    // dyn_hybrid_gr_init(G, src, dest, degree);
    // dyn_hybrid_gr_del(G, del_src, del_dest, degree, num_dels);
   
    /* Connectivity */

    /* Exp 1: Spanning forest construction */
    // dyn_gr_init(G, src, dest, degree);
    // dyn_st_forest_construct(G);
    
    /* Exp 2: Edge deletions and 
       spanning forest updates */
    // dyn_st_forest_ins(G, src, dest); 
    
    /* Exp 3: Connectivity Queries */
    // dyn_st_qeuries(G, del_src, del_dest, num_dels);

    /* Graph kernels */

    /* Exp 4: Induced Subgraphs */
    ts1 = 20;
    ts2 = 80;
    // dyn_induced_subgraphs(G, ts1, ts2);

    /* Exp 5: Path traversal */    
    
    // dyn_graph_traversal(G, src, dest);

    /* Exp 6: Path based centrality */

    
    free(filename);
    free(params);
    free(src);
    free(dest);
    free(del_src);
    free(del_dest);

    free(degree);
    free(G);
    free(taken);
    
    return 0;
}

void par_gen_RMAT_edges(dyn_graph_t* G, double* params, attr_id_t* src,
        attr_id_t* dest, attr_id_t* degree) {

    attr_id_t *permV;
    dyn_array_t* adj;
    attr_id_t *degree_hist;
    attr_id_t* edgeCount;
    attr_id_t *mem_chunk;
    long num_self_edges, num_dups;
    long num_mallocs, num_init_mallocs;
#ifdef _OPENMP
    omp_lock_t* vLock;
#endif
    attr_id_t edges_added, curr_num_edges_added;
    num_dups = num_self_edges = num_mallocs = num_init_mallocs = 0;

#ifdef _OPENMP
OMP("omp parallel reduction(+:num_dups) reduction(+:num_self_edges) \
    reduction(+:num_mallocs) reduction(+:num_init_mallocs) num_threads(50)")
#endif
{
    double el_time;
    long n, m;
    attr_id_t offset, avg_degree;
    long num_batches, num_edges_per_batch;
    
    long i, j;
    double a, b, c, d;
    double av, bv, cv, dv, S, p;
    double ins_del_ratio;
    int SCALE;
    double var;
    long step;
    long permute_vertices;
    attr_id_t tmpVal;
    attr_id_t* tmpArr;
    long u, v;
    int* stream, seed;
    attr_id_t brk;
    int tid, nthreads;
    attr_id_t *psrc, *pdest;
    attr_id_t pEdgeCount;
    attr_id_t num_edges_to_add;

#ifdef _OPENMP
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num(); 
#else
    nthreads = 1;
    tid = 0;
#endif

    if (tid == 0) 
        el_time = get_seconds();
   
    a = params[0];
    b = params[1];
    c = params[2];
    assert(a+b+c < 1);
    d = 1  - (a+b+c);
    permute_vertices = (long) params[3];
    num_batches = (long) params[4]; 
    ins_del_ratio = params[5];

    n = G->n;
    m = G->m;
    // num_edges_per_batch = m/num_batches;

    num_self_edges = 0;
    num_dups = 0;

    if (tid == 0) {
        mem_chunk = (attr_id_t *) malloc(2 * m  * sizeof(attr_id_t));
        adj = (dyn_array_t *) calloc(n, sizeof(dyn_array_t));
        degree_hist = (attr_id_t *) calloc(n, sizeof(attr_id_t));
#ifdef _OPENMP
        vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
#endif
        edgeCount = (attr_id_t *) malloc((nthreads+1)*sizeof(attr_id_t));
        
        assert(degree_hist != NULL);
        assert(mem_chunk != NULL);
        assert(adj != NULL);
#ifdef _OPENMP
        assert(vLock != NULL);
#endif
        assert(edgeCount != NULL);
    }
    
    avg_degree = (2*m)/n;

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
#ifdef _OPENMP
        omp_init_lock(&vLock[i]);
#endif
        adj[i].vals = &mem_chunk[i*avg_degree];
        adj[i].max_size = avg_degree;
    }

    if (tid == 0) {
        el_time = get_seconds() - el_time;
        fprintf(stderr, "Init time: %lf\n", el_time);
        el_time = get_seconds();
        fprintf(stderr, "Generating edges ... ");
        edges_added = 0; 
    }
 
    /* Initialize RNG stream */
    seed = 2387;
    stream = init_sprng(0, tid, nthreads, seed, SPRNG_DEFAULT);
    SCALE = log2(n);
    var = 0.1;

    num_edges_to_add = 1<<18;

    psrc = (attr_id_t *) malloc(num_edges_to_add * sizeof(attr_id_t));
    pdest = (attr_id_t *) malloc(num_edges_to_add * sizeof(attr_id_t));

#ifdef _OPENMP
OMP("omp barrier")
#endif

    /* Generate edges */
    while (edges_added < m) {

        if (num_edges_to_add > m - edges_added) {
            num_edges_to_add = m - edges_added;
        }

        pEdgeCount = 0;
#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for nowait")
#endif
        for (i=0; i<num_edges_to_add; i++) {

            u = 1;
            v = 1;
            step = n/2;

            av = a;
            bv = b;
            cv = c;
            dv = d;

            p = sprng(stream);
            if (p < av) {
                /* Do nothing */
            } else if ((p >= av) && (p < av+bv)) {
                v += step;
            } else if ((p >= av+bv) && (p < av+bv+cv)) {
                u += step;
            } else {
                u += step;
                v += step;
            }
        
            for (j=1; j<SCALE; j++) {
                step = step/2;

                /* Vary a,b,c,d by up to 10% */
                av *= 0.95 + var * sprng(stream);
                bv *= 0.95 + var * sprng(stream);
                cv *= 0.95 + var * sprng(stream);
                dv *= 0.95 + var * sprng(stream);

                S = av + bv + cv + dv;
                av = av/S;
                bv = bv/S;
                cv = cv/S;
                dv = dv/S;
            
                /* Choose partition */
                p = sprng(stream);
                if (p < av) {
                    /* Do nothing */
                } else if ((p >= av) && (p < av+bv)) {
                    v += step;
                } else if ((p >= av+bv) && (p < av+bv+cv)) {
                    u += step;
                } else {
                    u += step;
                    v += step;
                }
            }
       
            if (u == v) {
                num_self_edges++;
                continue;
            }

            u--;
            v--;
            brk = 0;
            if (adj[u].count < adj[v].count) {
#ifdef _OPENMP
                omp_set_lock(&vLock[u]);
#endif
                /* Check for duplicate in u's adjacencies */
                for (j = 0; j < adj[u].count; j++) {
                    if (adj[u].vals[j] == v) {
                        num_dups++;
                        break;
                    }
                }
                if (j < adj[u].count)
                    brk = 1;
#ifdef _OPENMP
                omp_unset_lock(&vLock[u]);
#endif
                if (brk)
                    continue;
           
            } else {
#ifdef _OPENMP
                omp_set_lock(&vLock[v]);
#endif
                /* Check for duplicate in v's adjacencies */
                for (j = 0; j < adj[v].count; j++) {
                    if (adj[v].vals[j] == u) {
                        num_dups++;
                        break;
                    }
                }
                if (j < adj[v].count)
                    brk = 1;
#ifdef _OPENMP
                omp_unset_lock(&vLock[v]);
#endif
                if (brk)
                    continue;
            }
        
            /* fprintf(stderr, "%ld %ld ", u, v); */
#ifdef _OPENMP
            omp_set_lock(&vLock[u]);
#endif
            if (adj[u].count == adj[u].max_size) {
                num_mallocs++;
                if (adj[u].max_size == avg_degree) {
                    num_init_mallocs++;
                    tmpArr = (attr_id_t *) malloc(4 * adj[u].max_size *
                         sizeof(attr_id_t));
                    assert(tmpArr != NULL);
                    memcpy(tmpArr, adj[u].vals, adj[u].max_size *
                        sizeof(attr_id_t));
                    adj[u].vals = NULL;
                    adj[u].vals = tmpArr;
                    adj[u].max_size = adj[u].max_size * 4;
                } else {
                    tmpArr = (attr_id_t *) malloc(4 * adj[u].max_size *
                         sizeof(attr_id_t));
                    assert(tmpArr != NULL);
                    memcpy(tmpArr, adj[u].vals, adj[u].max_size *
                        sizeof(attr_id_t));
                    free(adj[u].vals);
                    adj[u].max_size = adj[u].max_size * 4;
                    adj[u].vals = tmpArr;
                }
            }
            /* fprintf(stderr, "[u %ld %ld] ", adj[u].count, adj[u].max_size); */
            adj[u].vals[adj[u].count++] = v; 
#ifdef _OPENMP   
            omp_unset_lock(&vLock[u]);

            omp_set_lock(&vLock[v]);
#endif
            if (adj[v].count == adj[v].max_size) {
                num_mallocs++;
                if (adj[v].max_size == avg_degree) {
                    num_init_mallocs++;
                    tmpArr = (attr_id_t *) malloc(4 * adj[v].max_size *
                         sizeof(attr_id_t));
                    assert(tmpArr != NULL);
                    memcpy(tmpArr, adj[v].vals, adj[v].max_size *
                        sizeof(attr_id_t));
                    adj[v].vals = NULL;
                    adj[v].vals = tmpArr;
                    adj[v].max_size = adj[v].max_size * 4;
                } else {
                    tmpArr = (attr_id_t *) malloc(4 * adj[v].max_size *
                        sizeof(attr_id_t));
                    assert(tmpArr != NULL);
                    memcpy(tmpArr, adj[v].vals, adj[v].max_size *
                        sizeof(attr_id_t));
                    free(adj[v].vals);
                    adj[v].max_size = adj[v].max_size * 4;
                    adj[v].vals = tmpArr;
                }
            }
            /* fprintf(stderr, "[v %ld %ld] ", adj[v].count, adj[v].max_size); */
            adj[v].vals[adj[v].count++] = u;    
#ifdef _OPENMP
            omp_unset_lock(&vLock[v]);
#endif
            psrc[pEdgeCount] = u;
            pdest[pEdgeCount] = v;
            pEdgeCount++;
            // if (i % 10000 == 0) {
            //     fprintf(stderr, "%ld ", i);
            // }
        }

        edgeCount[tid+1] = pEdgeCount;
#ifdef _OPENMP
OMP("omp barrier")
#endif
        if (tid == 0) {
        
            edgeCount[0] = 0; 
            for (j=1; j<nthreads; j++) {
                edgeCount[j+1] += edgeCount[j];
            }
        }

#ifdef _OPENMP
OMP("omp barrier")
#endif
        /* Write partial results to src and dest */
        offset = edges_added + edgeCount[tid];
        for (j = 0; j<pEdgeCount; j++) {
            src[offset+j] = psrc[j];
            dest[offset+j] = pdest[j];
        }

#ifdef _OPENMP
OMP("omp barrier")
#endif
        if (tid == 0) {
            edges_added += edgeCount[nthreads];
            fprintf(stderr, "%ld ", edges_added);
        }
#ifdef _OPENMP
OMP("omp barrier")
#endif
    }

    free(psrc);
    free(pdest);

    if (tid == 0) {
        fprintf(stderr, "\n");
        el_time = get_seconds() - el_time;
        fprintf(stderr, "Edge generation time: %lf sec\n", el_time);
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (j = 0; j<n; j++) {
        if (adj[j].max_size != avg_degree)
            free(adj[j].vals);
    }

#ifdef _OPENMP
OMP("omp barrier")
#endif
    if (tid == 0) {
        free(mem_chunk);
        free(adj);
        free(edgeCount);
        el_time = get_seconds();
    }

#ifdef _OPENMP
OMP("omp barrier")
#endif

    if (permute_vertices) {

        if (tid == 0) {
            permV = (attr_id_t *) malloc(n*sizeof(attr_id_t));
            assert(permV != NULL);
        }
 
#ifdef _OPENMP
OMP("omp barrier")
#endif

#ifdef _OPENMP
OMP("omp for")
#endif       
        for (i=0; i<n; i++) {
            permV[i] = i;
        }

#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            j = n * sprng(stream);
            if (i != j) {
#ifdef _OPENMP
                omp_set_lock(&vLock[i]);
                omp_set_lock(&vLock[j]);
#endif
                tmpVal = permV[i];
                permV[i] = permV[j];
                permV[j] = tmpVal;
#ifdef _OPENMP
                omp_unset_lock(&vLock[j]);
                omp_unset_lock(&vLock[i]);
#endif
            }
        }

#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<m; i++) {
            // fprintf(stderr, "[%ld %ld] ", src[i], dest[i]);
            src[i]  = permV[src[i]];
            dest[i] = permV[dest[i]];
        }

        if (tid == 0) {
            el_time = get_seconds() - el_time;
            fprintf(stderr, "Permutation time: %lf sec\n", el_time);
            el_time = get_seconds();
        }
    }

#ifdef _OPENMP
OMP("omp for")
#endif
    for (i=0; i<m; i++) {
#ifdef _OPENMP
OMP("omp atomic")
#endif
        degree[src[i]]++;
    }
#if 0
#ifdef _OPENMP
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
#ifdef _OPENMP
OMP("omp atomic")
#endif
        degree_hist[degree[i]]++;
    }

    if (tid == 0) {
        for (i=0; i<n; i++) {
            if (degree_hist[i] != 0) {
                 fprintf(stderr, "%ld %ld\n", i, degree_hist[i]);
            }
        }
        el_time = get_seconds() - el_time;
        fprintf(stderr, "\nDegree hist time: %lf sec\n", el_time);
    }
#endif
#ifdef _OPENMP
OMP("omp barrier")
#endif

    if (tid == 0) {
#ifdef _OPENMP
        free(vLock);
#endif
        free(degree_hist);
        free(permV);
    }

#ifdef _OPENMP
OMP("omp barrier")
#endif

/*
OMP("omp for ")
    for (int i=0; i<n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
OMP("omp barrier")
*/
}

    fprintf(stderr, "num self edges: %ld, num dups: %ld\n", num_self_edges,
            num_dups);
    fprintf(stderr, "num mallocs: %ld, num init mallocs: %ld\n", num_mallocs,
            num_init_mallocs);
    
}

void read_dyn_test_config_file(dyn_graph_t* G, char* configfile, double* params) {

    /* read parameters from config file */
    FILE *fp;
    char line[128], var[32];
    double val;

    fp = fopen(configfile,"r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open config file: %s\n",configfile);
        exit(-1);
    }

    while (fgets(line, sizeof (line), fp) != NULL) {
        sscanf(line, "%s %lf", var, &val);
        if (*var == '#') continue;  /* comment */
        if (strcmp(var, "n") == 0) {
            G->n = (long) val;
            assert(G->n > 0);
        } else if (strcmp(var, "m") == 0) {
            G->m = (long) val;
            assert(G->m > 0);
        } else if (strcmp(var, "a") == 0) {
            params[0] = val;
            assert((val > 0) && (val < 1));
        } else if (strcmp(var, "b") == 0) {
            params[1] = val;
            assert((val > 0) && (val < 1));
        } else if (strcmp(var, "c") == 0) {
            params[2] = val;
            assert((val > 0) && (val < 1)); 
        } else if (strcmp(var, "permute_vertices") == 0) {
            params[3] = val;
        } else if (strcmp(var, "num_batches") == 0) {
            params[4] = val;
        } else if (strcmp(var, "ins_del_ratio") == 0) {
            params[5] = val;
        } else {
            fprintf(stderr,"Unknown parameter: %s\n", line);
        }
    }
    fclose(fp);

}

void dyn_ds_init_base(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    mem_chunk = (attr_id_t *) calloc(G->m, sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 11;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 4;
    nt[4] = 6; nt[5] = 8; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    long sum;
    double elt;
    int tid, nthreads;
    dyn_array_t* dt;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }
    sum = 0;
#ifdef _OPENMP
OMP("omp for")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        // mem_chunk[u]++;
        //dt = &G->adj[u];
        // dt->vals[dt->count++] = v;
        sum = sum + u + v;
        mem_chunk[u] = sum;
        // mem_chunk[u] = u + v;
        // degree_temp[u]++;
    }
    
    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        G->adj[i].count = 0;
    }
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(degree_temp);
    free(mem_chunk);
    free(G->adj);
#ifdef _OPENMP
    // free(vLock);
#endif

}

void dyn_ds_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 24;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 2;
    nt[4] = 4; nt[5] = 4; nt[6] = 8; nt[7] = 8;
    nt[8] = 12; nt[9] = 12; nt[9] = 16; nt[10] = 16;
    nt[11] = 24; nt[12] = 24; nt[14] = 32; nt[15] = 32;
    nt[16] = 40; nt[17] = 40; nt[18] = 48; nt[19] = 48;
    nt[20] = 56; nt[21] = 56; nt[22] = 64; nt[23] = 64;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for schedule(static)")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        
#ifdef _OPENMP        
        omp_set_lock(&vLock[u]);
        eid = G->adj[u].count++;
        omp_unset_lock(&vLock[u]);
#else
        eid = G->adj[u].count++;
#endif
        G->adj[u].vals[eid] = v;
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        G->adj[i].count = 0;
    }
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(degree_temp);
    free(mem_chunk);
    free(G->adj);
#ifdef _OPENMP
    // free(vLock);
#endif

}


void dyn_induced_subgraphs(dyn_graph_t* G, attr_id_t ts1, attr_id_t ts2) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;
    dyn_array_t* adj;
    attr_id_t* time_stamps;

    el_time = get_seconds();

}

void dyn_gr_del(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree, attr_id_t num_dels) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;
    dyn_array_t* adj;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    adj[0].vals = &mem_chunk[0];
    adj[0].max_size = degree[0];
    adj[0].count = degree[0];

    for (int i=1; i<G->n; i++) {
        adj[i].vals = &mem_chunk[degree_temp[i-1]];
        adj[i].max_size = degree[i];
        adj[i].count = degree[i];
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 11;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 4;
    nt[4] = 6; nt[5] = 8; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, j, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for schedule(static)")
#endif
    for (i = 0; i < num_dels; i++) {
        u = src[i];
        v = dest[i];
        
        for (j = 0; j<adj[u].max_size; j++) {
            if (adj[u].vals[j] == v) {
                adj[u].vals[j] = -1;        
#ifdef _OPENMP        
                omp_set_lock(&vLock[u]);
                adj[u].count--; 
                omp_unset_lock(&vLock[u]);
#else
                adj[u].count--;
#endif       
                break;
            }
        }
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, num_dels/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        adj[i].count = G->adj[i].count;
        for (j=0; j<adj[i].max_size; j++) {
            adj[i].vals[j] = G->adj[i].vals[j];
        }
    }

}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(degree_temp);
    free(mem_chunk);
    free(G->adj);
#ifdef _OPENMP
    // free(vLock);
#endif

}

void dyn_gr_init(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 1;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 24;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for schedule(static)")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        
#ifdef _OPENMP        
        omp_set_lock(&vLock[u]);
        eid = G->adj[u].count++;
        omp_unset_lock(&vLock[u]);
#else
        eid = G->adj[u].count++;
#endif
        G->adj[u].vals[eid] = v;
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(degree_temp);
#ifdef _OPENMP
    // free(vLock);
#endif

}

void dyn_st_forest_construct(dyn_graph_t* G) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    int chnk_sz;
    attr_id_t* S;
    long *start;
    char* visited;
    long phase_num, numPhases;
    long count;
    long *pSCount;
    attr_id_t *pred;
    attr_id_t src;

    el_time = get_seconds();
    
    S = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    visited = (char *) calloc(G->n, sizeof(char));
    start = (long *) calloc((numPhases+2), sizeof(long));
    pSCount = (long *) malloc(100*sizeof(long));
    pred = (attr_id_t *) malloc(G->n * sizeof(attr_id_t));

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    numPhases = 40;
 
    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);
    src = 0;

#ifdef _OPENMP
    nexp =10;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 4;
    nt[4] = 6; nt[5] = 8; nt[6] = 12; nt[7] = 16;
    nt[8] = 24; nt[9] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{

    attr_id_t *pS, *pSt;
    attr_id_t pCount, pS_size;
    attr_id_t v, w;
    int tid, nthreads;
    attr_id_t start_iter, end_iter;    
    attr_id_t j, k, vert, n;
    double elt;
    attr_id_t i;
#ifdef _OPENMP    
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
#else
    tid = 0;
    nthreads = 1;
#endif
    
    n = G->n;
    
    pS_size = (3*n)/nthreads + 1;
    pS = (attr_id_t *) malloc(pS_size*sizeof(attr_id_t));
    assert(pS != NULL);
#ifdef _OPENMP    
OMP("omp barrier")
#endif
    
#ifdef _OPENMP
OMP("omp barrier")
#endif
        
    if (tid == 0) {
        S[0] = src;
        visited[src] = (char) 1;
        count = 1;
        phase_num = 0;
        start[0] = 0;
        start[1] = 1;
    }

        
#ifdef _OPENMP
OMP("omp barrier")
#endif
    if (tid == 0) {
        elt = get_seconds();
    }    
    
    while (start[phase_num+1] - start[phase_num] > 0) {
        
            pCount = 0;

            start_iter = start[phase_num];
            end_iter = start[phase_num+1];
#ifdef _OPENMP
OMP("omp for nowait")
#endif
            for (vert=start_iter; vert<end_iter; vert++) {
            
                v = S[vert];
                for (j=0; j<G->adj[v].count; j++) {
                    w = G->adj[v].vals[j]; 
#ifdef _OPENMP
                    int myLock = omp_test_lock(&vLock[w]);
                    if (myLock) {
#endif
                        if (visited[w] != (char) 1) { 
                            visited[w] = (char) 1;
#if 0
                            if (pCount == pS_size) {
                                /* Resize pS */
                                pSt = (attr_id_t *)
                                    malloc(2*pS_size*sizeof(attr_id_t));
                                memcpy(pSt, pS, pS_size*sizeof(attr_id_t));
                                free(pS);
                                pS = pSt;
                                pS_size = 2*pS_size;
                            }
#endif
                            pS[pCount++] = w;
                            pred[w] = v;
                        }
#ifdef _OPENMP
                        omp_unset_lock(&vLock[w]);
                    }
#endif
                }
            }
            

            pSCount[tid+1] = pCount;
 
#ifdef _OPENMP
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
            
#ifdef _OPENMP
OMP("omp barrier")
#endif
            for (k = pSCount[tid]; k < pSCount[tid+1]; k++) {
                S[k] = pS[k-pSCount[tid]];
            } 
            
        } /* End of search */

    
        if (tid == 0) {
            elt = get_seconds() - elt;
            fprintf(stderr, "Number of threads: %d, Time taken: %lf seconds\n",
                    nthreads, elt);
            fprintf(stderr, "Search from vertex %ld, Number of vertices visited: %ld\n", src, count);
        }

#ifdef _OPENMP    
OMP("omp barrier")
    free(pS);
#endif
}
}


#ifdef _OPENMP
OMP("omp parallel for num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
 
    free(S);
    free(start);
    free(visited);
    free(pSCount);
#ifdef _OPENMP
    free(vLock);
#endif

    G->pred = pred;
#ifdef _OPENMP
    free(vLock);
#endif

}



void dyn_ds_init_malloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t avg_degree;
    attr_id_t num_mallocs, num_init_mallocs;
    attr_id_t* mem_chunk;
    int chnk_sz;

    el_time = get_seconds();

    mem_chunk = (attr_id_t *) malloc((2 * G->m + 1)  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    avg_degree = (2*G->m)/G->n;

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (long i=0; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[i*avg_degree];
        G->adj[i].max_size = avg_degree;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg 2 Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 6;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 32; nt[3] = 32;
    nt[4] = 64; nt[5] = 64; nt[6] = 8; nt[7] = 8;
    nt[8] = 12; nt[9] = 12; nt[9] = 16; nt[10] = 16;
    nt[11] = 24; nt[12] = 24; nt[14] = 32; nt[15] = 32;
    nt[16] = 40; nt[17] = 40; nt[18] = 48; nt[19] = 48;
    nt[20] = 56; nt[21] = 56; nt[22] = 64; nt[23] = 64;
    for (ti=0; ti<nexp; ti++) {
#ifdef MALLOC_COUNT
        num_mallocs = num_init_mallocs = 0;
#endif
OMP("omp parallel num_threads(nt[ti])") //reduction(+:num_mallocs) reduction(+:num_init_mallocs
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    attr_id_t* tmpArr;
    dyn_array_t* adju;

    double elt;
    int tid, nthreads;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 100;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif
    if (tid == 0)
        elt = get_seconds();

#ifdef _OPENMP
OMP("omp for schedule(dynamic, chnk_sz)")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
#ifdef _OPENMP
        omp_set_lock(&vLock[u]);
#endif
        adju = &G->adj[u];
        if (adju->count == adju->max_size) {
#ifdef MALLOC_COUNT
            num_mallocs++;
#endif
            if (adju->max_size == avg_degree) {
#ifdef MALLOC_COUNT
                num_init_mallocs++;
#endif
                tmpArr = (attr_id_t *) malloc(4 * adju->max_size *
                     sizeof(attr_id_t));
                // assert(tmpArr != NULL);
                memcpy(tmpArr, adju->vals, adju->max_size *
                    sizeof(attr_id_t));
                adju->vals = tmpArr;
                adju->max_size = adju->max_size * 4;
            } else {
                tmpArr = (attr_id_t *) malloc(4 * adju->max_size *
                     sizeof(attr_id_t));
                // assert(tmpArr != NULL);
                memcpy(tmpArr, adju->vals, adju->max_size *
                    sizeof(attr_id_t));
                free(adju->vals);
                adju->max_size = adju->max_size * 4;
                adju->vals = tmpArr;
            }
        }
        /* fprintf(stderr, "[u %ld %ld] ", adj[u].count, adj[u].max_size); */
        adju->vals[adju->count++] = v; 
#ifdef _OPENMP   
        omp_unset_lock(&vLock[u]);
#endif
   }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Number of threads: %ld, time: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        G->adj[i].count = 0;
        if (G->adj[i].max_size != avg_degree)
            free(G->adj[i].vals);
        G->adj[i].vals = &mem_chunk[i*avg_degree];
        G->adj[i].max_size = avg_degree;
    }
}

#ifdef MALLOC_COUNT
    fprintf(stderr, "num_mallocs: %ld, num_init_mallocs: %ld\n", num_mallocs,
            num_init_mallocs);
#endif
#ifdef _OPENMP
}
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif

    free(mem_chunk);
    free(G->adj);
#ifdef _OPENMP
    // free(vLock);
#endif

}

void dyn_vpart_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));

    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 6;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 32; nt[3] = 32;
    nt[4] = 64; nt[5] = 64; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;
    attr_id_t myt;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    elt = get_seconds();

    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        myt = u % nthreads;
        if (tid == myt) { 
            eid = G->adj[u].count++;
            G->adj[u].vals[eid] = v;
        }
    }

    elt = get_seconds() - elt;
    fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
            nthreads, elt, m/(elt*1000000));

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        G->adj[i].count = 0;
    }
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif

    free(mem_chunk);
    free(G->adj);
    free(degree_temp);
#ifdef _OPENMP
    // free(vLock);
#endif

}


void dyn_epart_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t* degree_temp;
    attr_id_t* degree_thresh;
    attr_id_t hd_count;
    attr_id_t hd_sum;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_thresh = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }
    
    hd_count = 0; 
    hd_sum = 0;

    for (int i=0; i<G->n; i++) {
        if (degree[i] > 1000) {
            degree_thresh[i] = 1;
            hd_count++;
            hd_sum += degree[i];
        }
    }

    fprintf(stderr, "Number of high degree vertices: %ld, adjacencies: %ld\n",
            hd_count, hd_sum);

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 6;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 32; nt[3] = 32;
    nt[4] = 64; nt[5] = 64; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    attr_id_t *adj1, *adj2;
    double elt;
    int tid, nthreads;
    attr_id_t count;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
#else
    tid = 0;
    nthreads = 1;
#endif

    adj1 = (attr_id_t *) malloc(((3*hd_sum)/nthreads)*sizeof(attr_id_t));
    adj2 = (attr_id_t *) malloc(((3*hd_sum)/nthreads)*sizeof(attr_id_t));
    count = 0;

    if (tid == 0) 
        elt = get_seconds();

#ifdef _OPENMP
OMP("omp for")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];
        if (degree_thresh[u]) {
            adj1[count] = u;
            adj2[count] = v;
            count++;
        } else {
#ifdef _OPENMP
            omp_set_lock(&vLock[u]);
            eid = G->adj[u].count++;
            omp_unset_lock(&vLock[u]);
#else
            eid = G->adj[u].count++;
#endif
            G->adj[u].vals[eid] = v;
        }
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
            nthreads, elt, m/(elt*1000000));
    }
    
#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        G->adj[i].count = 0;
    }

    free(adj1);
    free(adj2);
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif

    free(mem_chunk);
    free(G->adj);
    free(degree_temp);
    free(degree_thresh);
#ifdef _OPENMP
    // free(vLock);
#endif

}

void dyn_heap_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t *degree_thresh;
    attr_id_t hd_count;
    attr_id_t hd_sum;
    adj_bheap_t* Hadj;
    attr_id_t* Hval_mem_chunk;
        
    el_time = get_seconds();

    degree_thresh = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    hd_count = 0; 
    hd_sum = 0;

    for (int i=0; i<G->n; i++) {
        degree_thresh[i] = hd_sum;
        hd_count++;
        hd_sum += degree[i] + 1;
    }
    

    Hval_mem_chunk = (attr_id_t *) malloc((hd_sum + 1) * sizeof(attr_id_t));
    Hadj = (adj_bheap_t *) calloc(G->n, sizeof(adj_bheap_t));

    Hadj[0].vals = Hval_mem_chunk;
    for (int i=1; i<G->n; i++) {
        Hadj[i].vals = &Hval_mem_chunk[degree_thresh[i-1]];
    }
    
    fprintf(stderr, "Number of heaps created: %ld, memory: %ld\n",
            hd_count, hd_sum);

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 11;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 4;
    nt[4] = 6; nt[5] = 8; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;
    attr_id_t iv, jv, yv;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];

#ifdef _OPENMP      
            omp_set_lock(&vLock[u]);
#endif

            iv = ++(Hadj[u].count);
            while (iv >= 2) {
                jv = iv / 2;
                yv = Hadj[u].vals[jv];
                if (v >= yv) 
                    break;
                Hadj[u].vals[iv] = yv;
                iv = jv;
            }
            Hadj[u].vals[iv] = v;
#ifdef _OPENMP
            omp_unset_lock(&vLock[u]);
#endif
   }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        Hadj[u].count = 0;
    }
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(Hval_mem_chunk);
    free(Hadj);
    free(degree_thresh);
#ifdef _OPENMP
    // free(vLock);
#endif

}

void dyn_hybrid_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t *degree_temp, *degree_thresh;
    attr_id_t hd_count;
    attr_id_t hd_sum;
    adj_bheap_t* Hadj;
    attr_id_t* Hval_mem_chunk;
    attr_id_t first;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    el_time = get_seconds();

    degree_thresh = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    hd_count = 0; 
    hd_sum = 0;
    first = -1;
    
    for (int i=0; i<G->n; i++) {
        if (degree[i] > 64) {
            if (first == -1)
                 first = i;
            degree_thresh[i] = hd_sum;
            hd_sum += degree[i] + 1;
            hd_count++;
        }
    }
    

    Hval_mem_chunk = (attr_id_t *) malloc(2 * hd_sum * sizeof(attr_id_t));

    Hadj = (adj_bheap_t *) calloc(G->n, sizeof(adj_bheap_t));

    for (int i=0; i<G->n; i++) {
        if (degree_thresh[i]) {
            Hadj[i].vals = &Hval_mem_chunk[degree_thresh[i]];
        }
    }
    Hadj[first].vals = Hval_mem_chunk;

    degree_thresh[first] = 1;

    fprintf(stderr, "Number of high degree vertices: %ld, adjacencies: %ld\n",
            hd_count, hd_sum);

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 11;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 4;
    nt[4] = 6; nt[5] = 8; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;
    attr_id_t iv, jv, yv;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];

        if (degree_thresh[u]) {
#ifdef _OPENMP      
            omp_set_lock(&vLock[u]);
#endif

            iv = ++(Hadj[u].count);
            while (iv >= 2) {
                jv = iv / 2;
                yv = Hadj[u].vals[jv];
                if (v >= yv) 
                    break;
                Hadj[u].vals[iv] = yv;
                iv = jv;
            }
            Hadj[u].vals[iv] = v;
#ifdef _OPENMP
            omp_unset_lock(&vLock[u]);
#endif
        } else {        
#ifdef _OPENMP        
            omp_set_lock(&vLock[u]);
            eid = G->adj[u].count++;
            omp_unset_lock(&vLock[u]);
#else
            eid = G->adj[u].count++;
#endif
            G->adj[u].vals[eid] = v;
        }
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        G->adj[i].count = 0;
        Hadj[u].count = 0;
    }
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(degree_temp);
    free(degree_thresh);
    free(mem_chunk);
    free(Hval_mem_chunk);
    free(Hadj);
    free(G->adj);
#ifdef _OPENMP
    // free(vLock);
#endif

}


void dyn_hybrid_gr_init(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t *degree_temp, *degree_thresh;
    attr_id_t hd_count;
    attr_id_t hd_sum;
    adj_bheap_t* Hadj;
    attr_id_t* Hval_mem_chunk;
    attr_id_t first;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    el_time = get_seconds();

    degree_thresh = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    hd_count = 0; 
    hd_sum = 0;
    first = -1;

    for (int i=0; i<G->n; i++) {
        if (degree[i] > 64) {
            if (first == -1)
                first = i;
            degree_thresh[i] = hd_sum;
            hd_count++;
            hd_sum += degree[i] + 1;
        }
    }
    

    Hval_mem_chunk = (attr_id_t *) malloc((hd_sum + 1) * sizeof(attr_id_t));

    Hadj = (adj_bheap_t *) calloc(G->n, sizeof(adj_bheap_t));

    for (int i=0; i<G->n; i++) {
        if (degree_thresh[i]) {
            Hadj[i].vals = &Hval_mem_chunk[degree_thresh[i]];
        }
    }

    Hadj[first].vals = Hval_mem_chunk; 
    degree_thresh[first] = 1;
    
    fprintf(stderr, "Number of high degree vertices: %ld, adjacencies: %ld\n",
            hd_count, hd_sum);

    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    G->adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    G->adj[0].vals = &mem_chunk[0];
    G->adj[0].max_size = degree_temp[0];

    for (int i=1; i<G->n; i++) {
        G->adj[i].vals = &mem_chunk[degree_temp[i-1]];
        G->adj[i].max_size = degree_temp[i]-degree_temp[i-1];
        G->adj[i].count = 0;
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 1;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;
    attr_id_t iv, jv, yv;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for")
#endif
    for (i = 0; i < m; i++) {
        u = src[i];
        v = dest[i];

        if (degree_thresh[u]) {
#ifdef _OPENMP      
            omp_set_lock(&vLock[u]);
#endif

            iv = ++(Hadj[u].count);
            while (iv >= 2) {
                jv = iv / 2;
                yv = Hadj[u].vals[jv];
                if (v >= yv) 
                    break;
                Hadj[u].vals[iv] = yv;
                iv = jv;
            }
            Hadj[u].vals[iv] = v;
#ifdef _OPENMP
            omp_unset_lock(&vLock[u]);
#endif
        } else {        
#ifdef _OPENMP        
            omp_set_lock(&vLock[u]);
            eid = G->adj[u].count++;
            omp_unset_lock(&vLock[u]);
#else
            eid = G->adj[u].count++;
#endif
            G->adj[u].vals[eid] = v;
        }
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, m/(elt*1000000));
    }
    
    G->Hadj = Hadj;
    G->degree_thresh = degree_thresh;
}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif

    free(degree_temp);
    free(mem_chunk);
#ifdef _OPENMP
    // free(vLock);
#endif

}


void dyn_hybrid_gr_del(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree, attr_id_t num_dels) {

    double el_time;
#ifdef _OPENMP
    omp_lock_t* vLock;
    int nexp;
    int *nt, ti; 
#endif
    attr_id_t* mem_chunk;
    int chnk_sz;
    attr_id_t *degree_temp, *degree_thresh;
    attr_id_t hd_count;
    attr_id_t hd_sum;
    adj_bheap_t* Hadj;
    attr_id_t* Hval_mem_chunk;
    attr_id_t first;
    dyn_array_t* adj;

    el_time = get_seconds();

    degree_temp = (attr_id_t *) malloc(G->n*sizeof(attr_id_t));
    degree_thresh = (attr_id_t *) calloc(G->n, sizeof(attr_id_t));

    degree_temp[0] = degree[0];
    for (int i=1; i<G->n; i++) {
        degree_temp[i] = degree_temp[i-1] + degree[i];
    }

    hd_count = 0; 
    hd_sum = 0;
    first = -1;

    for (int i=0; i<G->n; i++) {
        if (degree[i] > 64) {
            if (first == -1)
                first = i;
            degree_thresh[i] = hd_sum;
            hd_count++;
            hd_sum += degree[i] + 1;
        }
    }

    Hval_mem_chunk = (attr_id_t *) malloc((hd_sum + 1) * sizeof(attr_id_t));

    Hadj = (adj_bheap_t *) calloc(G->n, sizeof(adj_bheap_t));

    for (int i=0; i<G->n; i++) {
        if (degree_thresh[i]) {
            Hadj[i].vals = &Hval_mem_chunk[degree_thresh[i]];
        }
    }
 
    Hadj[first].vals = Hval_mem_chunk; 
    degree_thresh[first] = 1;
    
    fprintf(stderr, "Number of high degree vertices: %ld, adjacencies: %ld\n",
            hd_count, hd_sum);


    mem_chunk = (attr_id_t *) malloc(G->m  * sizeof(attr_id_t));
    adj = (dyn_array_t *) calloc(G->n, sizeof(dyn_array_t));

    adj[0].vals = &mem_chunk[0];
    adj[0].max_size = degree[0];
    adj[0].count = degree[0];

    for (int i=1; i<G->n; i++) {
        adj[i].vals = &mem_chunk[degree_temp[i-1]];
        adj[i].max_size = degree[i];
        adj[i].count = degree[i];
    }

    el_time = get_seconds() - el_time;
    fprintf(stderr, "Alg Init time: %lf\n", el_time);

    el_time = get_seconds();
#ifdef _OPENMP
    vLock = (omp_lock_t *) malloc(G->n*sizeof(omp_lock_t));
#endif

#ifdef _OPENMP
OMP("omp parallel for schedule(static) num_threads(24)")
    for (long i=0; i<G->n; i++) {
        omp_init_lock(&vLock[i]);
    }
#endif
    el_time = get_seconds() - el_time;
    fprintf(stderr, "Lock init time: %lf\n", el_time);

#ifdef _OPENMP
    nexp = 11;
    nt = (int *) malloc(nexp*sizeof(int));
    nt[0] = 1; nt[1] = 1; nt[2] = 2; nt[3] = 4;
    nt[4] = 6; nt[5] = 8; nt[6] = 12; nt[7] = 16;
    nt[8] = 20; nt[9] = 24; nt[9] = 28; nt[10] = 32;
    for (ti=0; ti<nexp; ti++) {
OMP("omp parallel num_threads(nt[ti])")
#endif
{
    attr_id_t i, j, n, m;
    attr_id_t u, v, eid;
    double elt;
    int tid, nthreads;
    attr_id_t iv, jv, yv;
    attr_id_t nd;

    nd = num_dels;

    n = G->n;
    m = G->m;

#ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    chnk_sz = 1000;
OMP("omp barrier")
#else
    tid = 0;
    nthreads = 1;
#endif

    if (tid == 0) {
        elt = get_seconds();
    }

#ifdef _OPENMP
OMP("omp for")
#endif
    for (i = 0; i < nd; i++) {
        u = src[i];
        v = dest[i];

        if (degree_thresh[u]) {
#ifdef _OPENMP      
            omp_set_lock(&vLock[u]);
#endif

            iv = Hadj[u].count;
            iv--;
            while (iv >= 2) {
                jv = iv / 2;
                yv = Hadj[u].vals[jv];
                if (v >= yv) 
                    break;
                Hadj[u].vals[iv] = yv;
                iv = jv;
            }
            Hadj[u].vals[iv] = v;
#ifdef _OPENMP
            omp_unset_lock(&vLock[u]);
#endif
 

        } else {        
            for (j = 0; j<adj[u].max_size; j++) {
                if (adj[u].vals[j] == v) {
                    adj[u].vals[j] = -1;        
#ifdef _OPENMP        
                    omp_set_lock(&vLock[u]);
                    adj[u].count--; 
                    omp_unset_lock(&vLock[u]);
#else
                    adj[u].count--;
#endif       
                    break;
                }
            }
        }
    }

    if (tid == 0) {
        elt = get_seconds() - elt;
        fprintf(stderr, "Num threads: %d, time taken: %lf sec, TEPS: %lf\n",
                nthreads, elt, num_dels/(elt*1000000));
    }

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
    for (i=0; i<n; i++) {
        adj[i].count = G->adj[i].count;
        for (j=0; j<adj[i].max_size; j++) {
            adj[i].vals[j] = G->adj[i].vals[j];
        }
    }

}

#ifdef _OPENMP
   }
#endif

#ifdef _OPENMP
OMP("omp parallel for")
    for (long i=0; i<G->n; i++) {
        omp_destroy_lock(&vLock[i]);
    }
#endif
    free(degree_temp);
    free(mem_chunk);
    free(Hval_mem_chunk);
    free(Hadj);
    free(degree_thresh);
    free(G->adj);
#ifdef _OPENMP
    // free(vLock);
#endif

}


