#ifndef _GRAPH_GEN_H
#define _GRAPH_GEN_H
void graph_gen(graph_t* G, char* filename, char* graph_type);

void read_SNAP_graph(graph_t* G, char* filename);
void read_DIMACS_graph(graph_t* G, char* filename);
void read_METIS_graph(graph_t* G, char* filename);
void read_GraphML_graph(graph_t* G, char* filename);
void read_GML_graph(graph_t* G, char* filename);

void gen_RMAT_graph(graph_t* G, char* config_filename);
void read_RMAT_config_file(graph_t* G, char* configfile, double* params);

void read_dyn_test_config_file(dyn_graph_t* G, char* config_file, double*
        params);
void gen_random_graph(graph_t* G, char* config_filename);
void read_random_config_file(graph_t* G, char* configfile);
void gen_lmesh_graph(graph_t* G, char* config_filename);
void gen_sqmesh_graph(graph_t* G, char* config_filename);


#define GRAPH_GEN_SEED 34234 
void init_graph(graph_t* G);
void free_graph(graph_t* G);

#endif
