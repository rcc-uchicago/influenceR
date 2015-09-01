/* Header file for keyplayer.c */

#ifndef KEYPLAYER_H_
#define KEYPLAYER_H_

int keyplayer_driver(graph_t *g, int n, int k, double p, double tol, long maxsec, int *KP);
int keyplayer_driver_omp(graph_t *g, int n, int k, double p, double tol, long sec, long maxsec, int *KP);

#endif