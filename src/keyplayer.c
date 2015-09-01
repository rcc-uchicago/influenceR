/*
 keyplayer.c: Key Player implementation
 See also: http://www.bebr.ufl.edu/sites/default/files/Borgatti%20-%202006%20-%20Identifying%20sets%20of%20key%20players%20in%20a%20social%20networ.pdf
 
 AUTHOR: Simon Jacobs <sdjacobs@uchicago.edu>
 LICENSE: GPLv2
 
 Implementation of optimization part of Key Player algorithm. Includes a version
 parallelized with OpenMP.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <time.h>

#include "graph_defs.h"

#include "keyplayer-utils.h"

/* single core version: compute better solutions for T seconds */
/* Return code:
 * 0 if converges
 * 1 if does not converge 
*/
int keyplayer_driver(graph_t *g, int n, int k, double p, double tol, long maxsec, int *KP)
{
	
  int np, rank, new_rank = 0, stop;

  double *fits;
  int *allsets;
  time_t start;

  GetRNGstate(); // rather than srand
    
  problem_t problem;
  problem.graph = g;
  problem.distance = NULL;
  problem.round = 0;
  
  double fit, oldfit=-1;
  int ret = 1;

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
		
    if (fit - oldfit < tol) {
	    ret = 0;
		  break;
		}
		oldfit = fit;
	} while(difftime(time(0), start) < maxsec);
  
  for (int i = 0; i < n; i++)
    KP[i] = s[i];
  
  PutRNGstate();
  
	return ret;	
}

/* While we're working:
	 * compute better solutions for T seconds
	 * send my fit back to the master process.
	 * get a number back. if it's my rank, broadcast my s array to everyone! 
*/
int keyplayer_driver_omp(graph_t *g, int n, int k, double p, double tol, long maxsec, long sec, int *KP)
{
#ifndef _OPENMP
 int ret = keyplayer_driver(g, n, k, p, tol, maxsec, KP);
#else

  int np, rank, new_rank = 0, stop;

  double *fits;
  int *allsets;
  time_t start, fullstart;
  int ret = 1;

  
#pragma omp parallel shared(fits, allsets, new_rank, g, np, stop) private(rank, start, fullstart)
  {
  
    GetRNGstate(); // instead of srand
  
    
    np =  omp_get_num_threads();
    rank = omp_get_thread_num();
    
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
    double oldfit = -1;
	
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
  			
  			if (*fit - oldfit < tol) {
  			  ret = 0;
  			  break;
  			}
  			oldfit = *fit;
  			
  		} while(difftime(time(0), start) < sec);
      
		
  		//printf("Run %d, rank %d, fit %g\n", run, rank, *fit);
		
  		new_rank = 0;
		
#pragma omp barrier
      
  		/* Master process: find best fit. */
  		if (rank == 0) {
  			double max = 0;
  			for (int i = 0; i < np; i++) {
  				//printf("Run %d, rank %d, fit %g\n", run, i, fits[i]);
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
  
  PutRNGstate();
#endif
  
  return ret;
}
