#ifndef _GRAPH_KERNELS_H
#define _GRAPH_KERNELS_H
#include "graph_defs.h"

long graph_diameter(graph_t* G);
void strongly_connected_components(graph_t* G, attr_id_t* component_num);

attr_id_t biconnected_components(graph_t* G, attr_id_t* component_num);
void find_articulation_points(graph_t* G, attr_id_t* component_num);
void biconnected_components_recursive(graph_t* G, attr_id_t v, attr_id_t p_v);
void art_points_recursive(graph_t* G, attr_id_t u);
void BiCC_stack_push(attr_id_t a, attr_id_t b);
void BiCC_stack_pop(attr_id_t *a, attr_id_t *b);
void connected_components(graph_t* G, attr_id_t* component_num);
int aux_connected_components_init(graph_t* G);
int aux_connected_components_update(graph_t* G, attr_id_t mnum_components, attr_id_t maxbc_component);


long BFS_parallel_frontier_expansion(graph_t* G, long src, long diameter);
void BFS_path_limited_search(graph_t* G, long src, long diameter);
void BFS_sequential(graph_t* G, long src, int* d);
int vertex_cover_weighted(graph_t*);
int vertex_cover_unweighted(graph_t*);
#endif
