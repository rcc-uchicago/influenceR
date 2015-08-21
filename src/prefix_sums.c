/*
 AUTHOR: D.A. Bader and K. Madduri
 LICENSE: GPLv2
 See: http://snap-graph.sourceforge.net/

 Taken from utils.c
*/

#include "prefix_sums.h"

void prefix_sums(attr_id_t *input, attr_id_t* result, attr_id_t* p, attr_id_t n) {

#ifdef _OPENMP
    attr_id_t i, j, r, start, end, add_value;
    int tid, nthreads;

    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    r = n/nthreads;

    result[0] = 0;

OMP("omp for")
    for (i=1; i<n+1; i++)
        result[i] = input[i-1];

    start =  tid*r + 1;
    end   = (tid+1)*r;

    if (tid == nthreads-1)
        end = n+1;

    for (j=start; j<end; j++)
        result[j] = input[j-1] + result[j-1];

    p[tid] = result[end-1];

OMP("omp barrier")

    if (tid == 0) {
        for (j=1; j<nthreads; j++)
            p[j] += p[j-1];
    }

OMP("omp barrier")

    if (tid>0) {
        add_value=p[tid-1];
        for (j=start-1; j<end; j++)
            result[j] += add_value;
    }

OMP("omp barrier")
#else

    attr_id_t i;
    result[0] = 0;
    for (i=1; i<n+1; i++) {
        result[i] = result[i-1] + input[i-1];
    }
#endif

}
