/*
 keyplayer-utils.c: Utilities and functions for Key Player implementation
 See also: http://www.bebr.ufl.edu/sites/default/files/Borgatti%20-%202006%20-%20Identifying%20sets%20of%20key%20players%20in%20a%20social%20networ.pdf

 AUTHOR: Simon Jacobs <sdjacobs@uchicago.edu>
 LICENSE: GPLv2

 Contains a fast version of the KPP-NEG metric (metric 9 in Borgatti 2006).
*/

#include "graph_defs.h"
#include <Rmath.h>
#include <limits.h>

#include "keyplayer-utils.h"

#define X(i, j, n) ((i) * (n) + (j))

double * BFS_multiple(graph_t *g, int *src, int k, double *res);
double * BFS_single(graph_t *g, int src, double *res);
long BFS_parallel_frontier_expansion_with_distance(graph_t* G, long src, long diameter, double *distance);
void regen(int *gen, int *s, int *t, int n, int k);

int int_rand() {
  return (int) (unif_rand() * INT_MAX);
}

void gen_starting_set(int n, int k, int *s) {
    memset(s, 0, n * sizeof(int));
    for(int i = 0; i < k; i++) {
        int t = int_rand() % n;
        while(s[t] != 0)
            t = (t + 1) % n;
        s[t] = 1;
    }
}

/* D is a matrix where D[si,j] is the distance from s[si] to j. */
double kpmetric_graph(graph_t *g, double *D, int n, int *s, int *t, int k, int *argmin)
{
	double sum = 0;
	
	if (argmin != NULL)
		for (int i = 0; i < n; i++)
			argmin[i] = -1;
	
	for(int ti = 0; ti < n-k; ti++) {
		int i = t[ti];
		
		double min = INFINITY;
		
		for (int si = 0; si < k; si++) {
			double m = D[X(si, i, n)];
			if (m < min) {
				if (argmin != NULL)
					argmin[i] = si;
				min = m;
			}
		}
		
		if (min != 0 && min < INFINITY)
			sum += 1.0/min;
	}
	
	return sum/((double)n);
}
double kpmetric_graph_check(graph_t *g, double *D, int n, int *s, int *t, int k, int *prevargmin, int u, int v)
{
	
	double distance_v_to_all[n];
	BFS_single(g, v, distance_v_to_all);
	
	double sum = 0;
	
	for(int ti = 0; ti < n-k; ti++) {
		int i = t[ti];
		if (i == v)
			i = u;
		
		double min = INFINITY;
		
        int si = prevargmin[i];
        if (si != -1 && s[si] != u) {
            min = D[X(si, i, n)];
            double m = distance_v_to_all[i];
            if (m < min)
                min = m;
        }
        else {
            for (int si = 0; si < k; si++) {
                double m;
                if (s[si] != u)
                    m = D[X(si, i, n)];
                else
                    m = distance_v_to_all[i];
                if (m < min) {
                    min = m;
                }
            }
        }
		
		if (min != 0 && min < INFINITY)
			sum += 1.0/min;
	}
	
	return sum/((double)n);
}


/* Compute Metric 15 in Borghatti (19) */
/* D is a n-by-n matrix */
/* s is a list of indices, kp-set */
/* compute: X = sum(min(v not in s)(distance to w, w in s)) */
/* sum(1/X) / n */
double kpmetric_st(double *D, int n, int *s, int *t, int k, int *argmin)
{
	double sum = 0;

	if (argmin != NULL)
		for (int i = 0; i < n; i++)
			argmin[i] = -1;
	
    for(int ti = 0; ti < n-k; ti++) {
		
		int i = t[ti];
		
		double min = INFINITY;

		for (int si = 0; si < k; si++) {
			int j = s[si];
			double m = D[X(i, j, n)];
			if (m < min) {
				if (argmin != NULL)
					argmin[i] = j;
				min = m;
			}
		}
		
		if (min != 0 && min < INFINITY) /* This node is in s */
			sum += 1.0/((double)min);
    }
	
	return sum/((double) n);
}


/* Get next state function, probabilistically */
/* Method: pick a u randomly, and a v randomly. If fit improves, go for it. Otherwise, go for it with probability p. */
/* t = which(s) */
double get_next_state_graph(problem_t *this, int n, int *gen, int k, double p, int *ua, int *va, int round)
{
	int argmin[n];
	int s[k], t[n-k];
	regen(gen, s, t, n, k);
	
  graph_t *graph = this->graph;
  
  if (round != this->round) {
    free(this->distance);
    this->distance = NULL;
    this->round = round;
  }

	if (this->distance == NULL) {
		this->distance = (double *) malloc(n*k*sizeof(double));
		if (!this->distance) {
			return 0;
		}
		BFS_multiple(graph, s, k, this->distance);
	}

	double fit = kpmetric_graph(graph, this->distance, n, s, t, k, argmin);
	int tries = 0, u, v;
	while (1) {
		u = s[int_rand() % k];
#ifdef DEBUG
        // we want to match how we get a random v in the non-graph code. this sucks.
		v = int_rand % n;
		while(gen[v] != 0)
			v = (v + 1) % n;
#else
        v = t[int_rand() % (n-k)];
#endif				

		
  double fit_ = kpmetric_graph_check(graph, this->distance, n, s, t, k, argmin, u, v);
		
	if (fit_ > fit || unif_rand() < p) {
		*ua = u;
		*va = v;
		fit = fit_;
		break;
	}
	if (tries > 10000) {
		*ua = -1;
		*va = -1;
		break;
	}
	
	tries++;
}
	
	/* replace row u with row v */
	double * new_distance = (double *) malloc(n*k*sizeof(double));
	if (! new_distance) {
		//printf("error with malloc!");
		//exit(1);
    return 0;
	}
	int j = 0; // index into new_distance
	int vdone = 0;
	for (int i = 0; i < k; i++) {
		if (!vdone && s[i] > v) {
			BFS_single(graph, v, &new_distance[j*n]);
			j++;
			vdone = 1;
		}
		if (s[i] == u)
			continue;
		
		// copy
		for (int a = 0; a < n; a++)
			new_distance[X(j, a, n)] = this->distance[X(i, a, n)];
		j++;
	}
	if (j != k) {
		BFS_single(graph, v, &new_distance[j*n]);
		j++;
	}
	if (j != k) {
		//printf("fatal error!\n");
		//exit(1);
    return 0;
	}
		
	double *tmp = this->distance;
	this->distance = new_distance;
	free(tmp);
	
	return fit;
}




/* TODO: what's a good diameter value? */
double * BFS_multiple(graph_t *g, int *src, int k, double *res)
{
	int n = g->n;
	for (int i = 0; i < k*n; i++)
		res[i] = INFINITY;
	
	for (int i = 0; i < k; i++) {
		BFS_parallel_frontier_expansion_with_distance(g, (long) src[i], 75, &res[i*n]);
	}
	return res;
}

double * BFS_single(graph_t *g, int src, double *res)
{
	int n = g->n;
	for (int i = 0; i < n; i++)
		res[i] = INFINITY;
	BFS_parallel_frontier_expansion_with_distance(g, (long) src, 75, res);
	return res;
}

/* gen[i] = 0 if i in t, 1 if in s*/
void regen(int *gen, int *s, int *t, int n, int k)
{
                
	int si = 0, ti = 0;
	for (int i = 0; i < n; i++) {
		if(gen[i] == 1) {
			s[si] = i;
			si++;
		}
		else {
			t[ti] = i;  
			ti++;
		}
	}         
                        
	return;
}       


/* 
 Adapted from SNAP project, authors D.A. Bader and K. Madduri, see snap-graph.sourceforge.net
 See breadth_first_search.c : BFS_parallel_frontier_expansion
*/
long BFS_parallel_frontier_expansion_with_distance(graph_t* G, long src, long diameter, double *distance) {

    attr_id_t* S;
    long *start;
    char* visited;
    long *pSCount;
#ifdef DIAGNOSTIC
    double elapsed_time;
#endif
#ifdef _OPENMP_KP
    omp_lock_t* vLock;
#endif

    long phase_num, numPhases;
    long count;


#ifdef _OPENMP_KP

OMP("omp parallel")
    {
#endif

        attr_id_t *pS, *pSt;
        long pCount, pS_size;
        long v, w;
        int tid, nthreads;
        long start_iter, end_iter;    
        long j, k, vert, n;
#ifdef _OPENMP_KP
        int myLock;
#endif

#ifdef _OPENMP_KP    
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
#ifdef _OPENMP_KP
            vLock = (omp_lock_t *) malloc(n*sizeof(omp_lock_t));
#endif
        }

#ifdef _OPENMP_KP    
OMP("omp barrier")
OMP("omp for")
        for (i=0; i<n; i++) {
            omp_init_lock(&vLock[i]);
        }
#endif

#ifdef _OPENMP_KP
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


#ifdef _OPENMP_KP
OMP("omp barrier")
#endif

        while (start[phase_num+1] - start[phase_num] > 0) {

            pCount = 0;

            start_iter = start[phase_num];
            end_iter = start[phase_num+1];
#ifdef _OPENMP_KP
OMP("omp for")
#endif
            for (vert=start_iter; vert<end_iter; vert++) {

                v = S[vert];
                for (j=G->numEdges[v]; j<G->numEdges[v+1]; j++) {
                    w = G->endV[j]; 
                    if (v == w)
                        continue;
#ifdef _OPENMP_KP
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
#ifdef _OPENMP_KP
                        omp_unset_lock(&vLock[w]);
                    }
#endif
                }
            }


#ifdef _OPENMP_KP
OMP("omp barrier")
#endif            
            pSCount[tid+1] = pCount;

#ifdef _OPENMP_KP
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

#ifdef _OPENMP_KP
OMP("omp barrier")
#endif
            for (k = pSCount[tid]; k < pSCount[tid+1]; k++) {
                S[k] = pS[k-pSCount[tid]];
            } 


#ifdef _OPENMP_KP
OMP("omp barrier")
#endif
        } /* End of search */

#ifdef DIAGNOSTIC
        if (tid == 0) {
            REprintf("Search from vertex %ld," 
                    " No. of vertices visited: %ld\n", src, count);
        }
#endif

        free(pS);
#ifdef _OPENMP_KP    
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
#ifdef _OPENMP_KP
            free(vLock);
#endif

        }

#ifdef _OPENMP_KP    
    }
#endif

#ifdef DIAGNOSTIC    
    elapsed_time = get_seconds() - elapsed_time;
    REprintf("Time taken for BFS: %lf seconds\n", elpased_time);
#endif
    return count;
}
