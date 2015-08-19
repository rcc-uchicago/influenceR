#ifndef _GRAPH_METRICS_H
#define _GRAPH_METRICS_H

void vertex_betweenness_centrality(graph_t* G, double* BC, long numSrcs);
void vertex_betweenness_centrality_simple(graph_t* G, double* BC, long numSrcs);
void vertex_betweenness_centrality_parBFS(graph_t* G, double* BC, long numSrcs);


#endif
