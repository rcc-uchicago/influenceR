#ifndef _UTILS_H
#define _UTILS_H
#include "graph_defs.h"

/* utils.c */
double get_seconds(void);
void prefix_sums(attr_id_t*, attr_id_t*, attr_id_t*, attr_id_t);
void usage_graph_options(void);
void graph_ext_check(char*, char*);
void print_graph(graph_t*);
void print_snap_header(FILE*);
void print_graph_header(FILE*, graph_t*, const char*);
void save_undir_unwgt_graph(FILE*, graph_t*);
void generate_random_walk_seeds(graph_t*, int, int, attr_id_t*);

/* vectorUtils.c */

void printDoubleVector(double*,int,int);
void printIntVector(int*,int,int);
void print_attr_id_t_Vector(attr_id_t*,attr_id_t,attr_id_t);

/* list.c */
list_t* makeList(void);
node_t* makeNode(int);
void append(list_t*,node_t*);
node_t *getFirst(list_t*);
void deleteFirst(list_t*);
void printList(list_t*);


#endif
