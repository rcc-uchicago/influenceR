/* Header file for keyplayer-utils.c */

#ifndef KEYPLAYER_UTILS_H
#define KEYPLAYER_UTILS_H

typedef struct problem_struct {
  graph_t *graph;
  int round;
  double *distance;
} problem_t;

void gen_starting_set(int n, int k, int *s);
double get_next_state_graph(problem_t *this, int n, int *gen, int k, double p, int *ua, int *va, int round);

#endif
