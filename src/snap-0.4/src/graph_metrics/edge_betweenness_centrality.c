#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"
#include "sprng.h"

#define CACHELOG 7
#define NOSHARE(x) ((x)<<CACHELOG)

void evaluate_edge_centrality_bcpart(graph_t* g, attr_id_t* ebc_eval_data1, 
        double* ebc_eval_data2, comm_list_bc_t* comm_list, 
        attr_id_t num_components, 
        attr_id_t curr_component1, attr_id_t curr_component2) {


    /* stack of vertices in the order of non-decreasing 
     * distance from s. Also used to implicitly 
     * represent the BFS queue. 
     */
    attr_id_t* visited_vertices;    

    /* A struct to store information about the visited vertex */
    typedef struct {
#ifdef _OPENMP
        omp_lock_t vlock;
#endif
        attr_id_t sigma;         /* No. of shortest paths */
        attr_id_t child_count;   /* Child count associated with each vertex */ 
        double delta;            /* Partial dependency    */
        unsigned short d;        /* distance of vertex from source vertex */		
        /* char padding[32];*/	 /* If necessary, we can pad this structure to
                                    align to cache line size */
    } S_info_t;

    S_info_t *vis_path_info;

    attr_id_t* child;      /* vertices that succeed v on 
                              shortest paths from s */
    attr_id_t* child_epos; /* Edge corresponding to this child vertex */
    attr_id_t* vis_srcs;   /* Randomly chosen set of vertices 
                              to initiate traversals from */
#ifdef _OPENMP
    omp_lock_t* vLock;     /* Need locks for permuting vertices */
#endif

    long* num_visited;
    long* visited_counts;
    double elapsed_time;
#ifdef _OPENMP    
OMP("omp parallel")
#endif
    {
        /* Vars to index local thread stacks */
        attr_id_t* p_vis;
        attr_id_t* p_vis_temp;
        long p_vis_size, p_vis_count;
        long phase_num;

        /* Other vars */
        long i, j, k, p, s;
        long v, w;

        long n, m, num_traversals;
        c_vert_t* cvl;
        c_edge_t* cel;
        long d_phase;
        int tid, nthreads;
        int* stream;
        long MAX_NUM_PHASES;
        int SPRNG_SEED;
        long ph_start, ph_end;
        long sigma_v;
        attr_id_t num_edges_v;
        long offset;

        int seed;

        S_info_t *vis_path_info_i, *vis_path_info_v, *vis_path_info_w;
        attr_id_t vis_path_info_v_child_count;
        attr_id_t vis_path_info_v_sigma;
        double vis_path_info_v_delta;
        int d_val;
        double incr;
        attr_id_t epos;

#if DIAGNOSTIC
        double elapsed_time_part, elapsed_time_all;
#endif

#ifdef _OPENMP
        int myLock, l1, l2;
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif

        SPRNG_SEED = 129832;

#if DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_all = get_seconds();
            elapsed_time_part = get_seconds();
        }
#endif

        n = g->n;
        m = g->m;
        cvl = g->cvl;
        cel = g->cel;
        /* Generate a random stream of source vertices */
        if (tid == 0) {
            vis_srcs = (attr_id_t *) malloc(n*sizeof(attr_id_t));
            assert(vis_srcs != NULL);
#ifdef _OPENMP
            vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
            assert(vLock != NULL);
#endif
        }

        /* Initialize RNG stream */
        seed = SPRNG_SEED;
        stream = init_sprng(SPRNG_CMRG, tid, nthreads, seed, SPRNG_DEFAULT);

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            vis_srcs[i] = i;
#ifdef _OPENMP
            omp_init_lock(&vLock[i]);
#endif
        }

#ifdef _OPENMP
OMP("omp for")
#endif    
        for (i=0; i<n; i++) {
            j = n*sprng(stream);
            if (i != j) {
#ifdef _OPENMP
                l1 = omp_test_lock(&vLock[i]);
                if (l1) {
                    l2 = omp_test_lock(&vLock[j]);
                    if (l2) {
#endif
                        k = vis_srcs[i];
                        vis_srcs[i] = vis_srcs[j];
                        vis_srcs[j] = k;
#ifdef _OPENMP
                        omp_unset_lock(&vLock[j]);
                    }
                    omp_unset_lock(&vLock[i]);
                }
#endif
            }
        }

        /* Clear current betweenness values */
        if (curr_component1 != -1) {
#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for schedule(guided)")
#endif
            for (i=0; i<n; i++) {
                if ((cvl[i].comm_id != curr_component1) && 
                        (cvl[i].comm_id != curr_component2))
                    continue;
                for (j=cvl[i].num_edges; j<cvl[i+1].num_edges; j++) {
                    cel[j].cval = 0.0;
                }

            }
        }

        if (tid == 0) {
#ifdef _OPENMP
            free(vLock);
#endif
#if DIAGNOSTIC
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stdout, "Vertex ID permutation time: %lf seconds\n", 
                    elapsed_time_part);
            elapsed_time_part = get_seconds();
#endif
        }
        /* Start timing code from here */

        if (tid == 0) {
            elapsed_time = get_seconds();
        }

        MAX_NUM_PHASES = 1000;

        if (tid == 0) {
            /* Allocate memory for the data structures */
            visited_vertices = (attr_id_t *) malloc(n * sizeof(attr_id_t));
            vis_path_info = (S_info_t *) malloc(n * sizeof(S_info_t));
            child = (attr_id_t *) malloc(m * sizeof(attr_id_t));
            child_epos = (attr_id_t *) malloc(m * sizeof(attr_id_t));
            assert(visited_vertices != NULL);
            assert(vis_path_info != NULL);
            assert(child != NULL);
            assert(child_epos != NULL);
            num_visited = (long *) malloc(MAX_NUM_PHASES*sizeof(long));
            visited_counts = (long *) malloc(NOSHARE(nthreads+1)*sizeof(long));

            assert(num_visited != NULL);
            assert(visited_counts != NULL);
        }

        /* local memory for each thread */
        p_vis_size = (2*n)/nthreads;
        p_vis = (attr_id_t *) malloc(p_vis_size*sizeof(attr_id_t));
        assert(p_vis != NULL);
        num_traversals = 0;
        p_vis_count = 0;

#if DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stdout, "Memory allocation time: %9.6lf seconds\n", 
                    elapsed_time_part);
            elapsed_time_part = get_seconds();
        }
#endif

#ifdef _OPENMP
OMP("omp barrier")
#endif

#ifdef _OPENMP
OMP("omp for")
#endif
        for (i=0; i<n; i++) {
            vis_path_info_i = &(vis_path_info[i]);
            vis_path_info_i->d = 0;
            vis_path_info_i->sigma = 0;
            vis_path_info_i->child_count = 0;
#ifdef _OPENMP
            omp_init_lock(&(vis_path_info_i->vlock));
#endif
        }

#if DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_part = get_seconds() - elapsed_time_part;
            fprintf(stdout, "BC initialization time: %9.6lf seconds\n", 
                    elapsed_time_part);
        }
#endif

        /* Start traversals */
        for (p = 0; p < n; p ++) {

            s = vis_srcs[p];
            if ((cvl[s+1].num_edges - cvl[s].num_edges) == 0) {
                continue;
            }

            if (curr_component1 != -1) {
                if ((cvl[s].comm_id != curr_component1) && 
                        (cvl[s].comm_id != curr_component2))
                    continue;
            }

            num_traversals++;

            /*
               if (num_traversals == numV + 1) {
               break;
               }
             */

            phase_num = 0;

            if (tid == 0) {

                visited_vertices[0] = s;
                vis_path_info[s].sigma = 1;
                vis_path_info[s].d = 1;
                vis_path_info[s].child_count = 0;

                num_visited[phase_num] = 0;
                num_visited[phase_num+1] = 1;
                visited_counts[NOSHARE(0)] = 0;
            }

            visited_counts[NOSHARE(tid+1)] = 0;

            d_phase = 1;

#if DIAGNOSTIC
            if (tid == 0) 
                elapsed_time_part = get_seconds();
#endif

#ifdef _OPENMP       
OMP("omp barrier")
#endif

            while (num_visited[phase_num+1] - num_visited[phase_num] > 0) {

                ph_start = num_visited[phase_num];
                ph_end = num_visited[phase_num+1];

                p_vis_count = 0;
                d_phase++;

#ifdef _OPENMP
OMP("omp barrier")
OMP("omp for schedule(guided,4)")
#endif
                for (i = ph_start; i < ph_end; i++) {
                    v = visited_vertices[i];
                    vis_path_info_v = &(vis_path_info[v]);
                    sigma_v = vis_path_info_v->sigma;
                    num_edges_v = cvl[v].num_edges;

                    for (j = num_edges_v; j < cvl[v+1].num_edges; j++) {

                        /* Filter edges that are deleted */
                        if (cel[j].mask == 0) {

                            w = cel[j].dest;
                            vis_path_info_w = &(vis_path_info[w]);

                            d_val = vis_path_info_w->d;
                            if (d_val == 0) {
#ifdef _OPENMP
                                myLock = omp_test_lock(&(vis_path_info_w->vlock));
                                if (myLock == 1) {
#endif
                                    vis_path_info_w->d = d_phase;
                                    vis_path_info_w->sigma += sigma_v;
#ifdef _OPENMP
                                    omp_unset_lock(&(vis_path_info_w->vlock));
#endif
                                    child[num_edges_v+vis_path_info_v->child_count] = w;
                                    child_epos[num_edges_v + 
                                        vis_path_info_v->child_count] = j;
                                    vis_path_info_v->child_count++;
                                    /* Add w to local stack */
                                    if (p_vis_count == p_vis_size) {
                                        /* Resize p_vis */
                                        p_vis_temp = (attr_id_t *)
                                            malloc(2*p_vis_size*sizeof(attr_id_t));
                                        assert(p_vis_temp != NULL);
                                        memcpy(p_vis_temp, p_vis, 
                                                p_vis_size*sizeof(attr_id_t));
                                        free(p_vis);
                                        p_vis = p_vis_temp;
                                        p_vis_size = 2*p_vis_size;
                                    }
                                    p_vis[p_vis_count++] = w;
#ifdef _OPENMP
                                } else {
                                    child[num_edges_v + 
                                        vis_path_info_v->child_count] = w;
                                    child_epos[num_edges_v + 
                                        vis_path_info_v->child_count] = j;
                                    vis_path_info_v->child_count++;

                                    omp_set_lock(&(vis_path_info_w->vlock));
                                    vis_path_info_w->sigma += sigma_v;
                                    omp_unset_lock(&(vis_path_info_w->vlock));
                                }	
#endif

                            } else if (d_val == d_phase) {
#ifdef _OPENMP
                                omp_set_lock(&(vis_path_info_w->vlock));
#endif
                                vis_path_info_w->sigma += sigma_v;
#ifdef _OPENMP
                                omp_unset_lock(&(vis_path_info_w->vlock));
#endif
                                child[num_edges_v + 
                                    vis_path_info_v->child_count] = w;
                                child_epos[num_edges_v + 
                                    vis_path_info_v->child_count] = j;
                                vis_path_info_v->child_count++;

                            }

                        } /* End of edge filter loop */

                    } /* End of inner loop */
                } /* End of outer for loop */

                /* Merge all local stacks for next iteration */
                phase_num++;

                visited_counts[NOSHARE(tid+1)] = p_vis_count;
                /* fprintf(stdout, "s: %d, tid: %d, phase %d, count %d\n", 
                   s, tid,
                 * phase_num, p_vis_count); */
#ifdef _OPENMP
OMP("omp barrier")
#endif

                if (tid == 0) {
                    visited_counts[NOSHARE(0)] = num_visited[phase_num];
                    for(k=1; k<=nthreads; k++) {
                        visited_counts[NOSHARE(k)] += visited_counts[NOSHARE(k-1)];
                    }
                    num_visited[phase_num+1] = visited_counts[NOSHARE(nthreads)];
                }


#ifdef _OPENMP           
OMP("omp barrier")
#endif

                for (k = visited_counts[NOSHARE(tid)]; 
                        k < visited_counts[NOSHARE(tid+1)]; k++) {
                    offset = visited_counts[NOSHARE(tid)];
                    visited_vertices[k] = p_vis[k-offset];
                } 


#ifdef _OPENMP            
OMP("omp barrier")
#endif
            }

#if DIAGNOSTIC 
            if (tid == 0) {
                elapsed_time_part = get_seconds() - elapsed_time_part;
                fprintf(stdout, "Traversal time: %9.6lf seconds, src"
                        " %d, num phases %d\n",elapsed_time_part, s, phase_num);
            }
#endif

            phase_num--;

            /*        
             if (tid == 0)  
                fprintf(stdout, "s: %d, phase %d, start %d, end %d\n", 
                        s, phase_num,
                        num_visited[phase_num], num_visited[phase_num+1]);
            */
            /* Clear delta scores of leaf vertices */
            if (tid == 0) {
                for (j = num_visited[phase_num]; 
                        j < num_visited[phase_num+1]; j++) {
                    v = visited_vertices[j];
                    vis_path_info[v].delta = 0.0;
                }
            }

            phase_num--;

#if DIAGNOSTIC
            if (tid ==0)
                elapsed_time_part = get_seconds();
#endif

#ifdef _OPENMP        
OMP("omp barrier")
#endif
            /* Dependency accumulation phase */
            while (phase_num >= 0) {
                ph_start = num_visited[phase_num];
                ph_end = num_visited[phase_num+1];
#ifdef _OPENMP
OMP("omp for schedule(guided,4) nowait")
#endif
                for (j = ph_start; j < ph_end; j++) {
                    v = visited_vertices[j];
                    vis_path_info_v = &vis_path_info[v];
                    vis_path_info_v_child_count = vis_path_info_v->child_count;
                    vis_path_info_v_delta = 0.0;
                    vis_path_info_v_sigma = vis_path_info_v->sigma;
                    num_edges_v = cvl[v].num_edges;
                    for (k = 0; k < vis_path_info_v_child_count; k++) {
                        w = child[num_edges_v+k];
                        epos = child_epos[num_edges_v+k];
                        incr = vis_path_info_v_sigma * 
                            (1+vis_path_info[w].delta)/vis_path_info[w].sigma;
                        cel[epos].cval += incr;
                        vis_path_info_v_delta += incr;
                    }
                    vis_path_info_v->delta = vis_path_info_v_delta;
                }

                phase_num--;

#ifdef _OPENMP
OMP("omp barrier")
#endif            
            }

#if DIAGNOSTIC
            if (tid == 0) {
                elapsed_time_part = get_seconds() - elapsed_time_part;
                fprintf(stdout, "Accumulation time: %9.6lf seconds\n", 
                        elapsed_time_part);
            }
#endif

#ifdef _OPENMP
OMP("omp for nowait")
#endif
            for (j=0; j<n; j++) {
                vis_path_info[j].d = 0;
                vis_path_info[j].sigma = 0;
                vis_path_info[j].child_count = 0;
            }

#ifdef _OPENMP
OMP("omp barrier")
#endif
        }

#if DIAGNOSTIC
        if (tid == 0) {
            elapsed_time_all = get_seconds() - elapsed_time_all;
            fprintf(stdout, "BC computation time: %9.6lf seconds\n",
                    elapsed_time_all);
        }
#endif

        free(p_vis);

        if (tid == 0) { 
            free(visited_vertices);
            free(vis_path_info);
            free(num_visited);
            free(visited_counts);
            free(child);
            free(child_epos);
            elapsed_time = get_seconds() - elapsed_time;
            free(vis_srcs);
        }

        free_sprng(stream);

#ifdef _OPENMP
OMP("omp barrier")
#endif
    }

}    

