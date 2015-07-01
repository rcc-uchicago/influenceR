#ifndef KEYPLAYER_H_
#define KEYPLAYER_H_

void keyplayer_driver(graph_t *g, int n, int k, double p, double tol, long maxsec);
void keyplayer_driver_parallel(graph_t *g, int n, int k, double p, double tol, long sec, long maxsec);

#endif