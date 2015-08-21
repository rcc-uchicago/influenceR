/* Header file for bridging.c */

#include "graph_defs.h"

#ifndef BRIDGING_H
#define BRIDGING_H

double *bridging(graph_t *G, int *edgelist, double *scores);

#ifdef USE_MPI
double *bridging_MPI(graph_t *G, int *edgelist, double *scores);
#endif

#endif