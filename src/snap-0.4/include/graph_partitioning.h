#if !defined(GRAPH_PARTITIONING_H_)
#define GRAPH_PARTITIONING_H_
#include "graph_defs.h"
void modularity_spectral(graph_t *G, attr_id_t *membership, 
        attr_id_t *numCommunities, attr_id_t use_improvement);
void modularity_spectral_wo_klin(graph_t *G, attr_id_t *membership, 
        attr_id_t *numCommunities);
void computeModularityValue(graph_t *G, attr_id_t *membership, 
        attr_id_t numCommunities, double *modularity);
void modularity_betweenness(graph_t *g, attr_id_t *membership, 
        attr_id_t *num_communities, double *modularity, 
        double sampling_val);
void modularity_greedy_agglomerative(graph_t *g, char *alg_type, 
        attr_id_t *membership, attr_id_t *num_communities, double *modularity);
void seed_set_community_detection(graph_t*, char*, attr_id_t*, attr_id_t,
				  attr_id_t*, double*);
void pagerank_community (graph_t*, const attr_id_t*, attr_id_t,
			 attr_id_t*, double, double);
void andersen_lang (graph_t*, const attr_id_t*, attr_id_t, attr_id_t*, const int);
void BFS_seed_set_expansion (graph_t*, attr_id_t*, int, attr_id_t*, attr_id_t);
#endif /* GRAPH_PARTITIONING_H_ */
