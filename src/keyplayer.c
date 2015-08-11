#include <R.h>
#include <Rinternals.h>

#include <time.h>

#include <graph_defs.h>

#include "keyplayer-utils.h"

/* single core version: compute better solutions for T seconds */
void keyplayer_driver(graph_t *g, int n, int k, double p, double tol, long maxsec, int *KP)
{
	
  int np, rank, new_rank = 0, stop;

  double *fits;
  int *allsets;
  time_t start;

  srand(time(NULL));
    
  problem_t problem;
  problem.graph = g;
  problem.distance = NULL;
  problem.round = 0;
  
  double fit;

	int s[n];
	gen_starting_set(n, k, s);

  start = time(0);
	do {
		int u, v;

		fit = get_next_state_graph(&problem, n, s, k, p, &u, &v, 0);
		
		if (u >= 0)
			s[u] = 0;
		if (v >= 0)
			s[v] = 1;
	} while(difftime(time(0), start) < maxsec);
  
  for (int i = 0; i < n; i++)
    KP[i] = s[i];
  
	return;	
}

/* While we're working:
	 * compute better solutions for T seconds
	 * send my fit back to the master process.
	 * get a number back. if it's my rank, broadcast my s array to everyone! 
*/
void keyplayer_driver_omp(graph_t *g, int n, int k, double p, double tol, long maxsec, long sec, int *KP)
{
#ifndef OPENMP
 keyplayer_driver(g, n, k, p, tol, maxsec, KP);
#else
  int np, rank, new_rank = 0, stop;

  double *fits;
  int *allsets;
  time_t start, fullstart;

#pragma omp parallel shared(fits, allsets, new_rank, g, np, stop) private(rank, start, fullstart)
  {
    
    np =  omp_get_num_threads();
    rank = omp_get_thread_num();
    
    srand(time(NULL) + rank);
    
    if (rank == 0) {
      allsets = (int *) R_alloc(n * np, sizeof(int));
      fits = (double *) R_alloc(np, sizeof(double));
    }
    
    problem_t problem;
    problem.graph = g;
    problem.distance = NULL;
    problem.round = 0;
    
#pragma omp barrier
    
  	int *s = &allsets[rank * n];
  	gen_starting_set(n, k, s);
	
  	start = time(0), fullstart = start;
  	double *fit = &fits[rank];
    *fit = 0;
    double oldfit = 0;
	
  	int run = 0;
	
  	do {
		
  		start = time(0);
		
  		do {
  			int u, v;
		
  			*fit = get_next_state_graph(&problem, n, s, k, p, &u, &v, run);
				
  			if (u >= 0)
  				s[u] = 0;
  			if (v >= 0)
  				s[v] = 1;
  		} while(difftime(time(0), start) < sec);
      
		
  		//printf("Run %d, rank %d, fit %g\n", run, rank, *fit);
		
  		new_rank = 0;
		
#pragma omp barrier
      
  		/* Master process: find best fit. */
  		if (rank == 0) {
  			double max = 0;
  			for (int i = 0; i < np; i++) {
  				printf("Run %d, rank %d, fit %g\n", run, i, fits[i]);
  				if (fits[i] > max) {
  					max = fits[i];
  					new_rank = i;
  				}
  			}
  			if (max - oldfit < tol || (difftime(time(0), fullstart) > maxsec)) {
  				stop = 1;
  			}
  			oldfit = max;
  		}

#pragma omp barrier
		
  		/* update s, or send it */
      if (rank != new_rank) {
        int *best_s = &allsets[n * new_rank];
        for (int i = 0; i < n; i++)
          s[i] = best_s[i];
      }		
		
  	  run++;
      
#pragma omp barrier
      
  	} while(!stop);
  }
  
  int *s = &allsets[0]; // s set for rank 0, which contains the best (as do all &allsets[i*n])
  for (int i = 0; i < n; i++)
    KP[i] = s[i];
#endif
}
