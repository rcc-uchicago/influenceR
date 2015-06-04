#include "graph_defs.h"
#include "graph_gen.h"

void graph_gen(graph_t* g, char* filename, char* graph_type) {

    if (strcmp(graph_type, "snap") == 0) {
        read_SNAP_graph(g, filename);
    } else if (strcmp(graph_type, "dimacs") == 0) {
        read_DIMACS_graph(g, filename);
    } else if (strcmp(graph_type, "metis") == 0) {
        read_METIS_graph(g, filename);
    } else if (strcmp(graph_type, "gml") == 0) {
        read_GML_graph(g, filename);
    } else if (strcmp(graph_type, "rand") == 0) {
        gen_random_graph(g, filename);
    } else if (strcmp(graph_type, "rmat") == 0) {
        gen_RMAT_graph(g, filename);
    } else if (strcmp(graph_type, "sqm") == 0) {
        gen_sqmesh_graph(g, filename);
    } else if (strcmp(graph_type, "lm") == 0) {
        gen_lmesh_graph(g, filename);
    } else {
        fprintf(stderr, "Error! Invalid graph format (%s). Exiting ...\n", 
                graph_type);
        exit(-1);
    }

}

void init_graph(graph_t* g) {

    g->n = 0;
    g->m = 0;
    g->weight_type = -1;
    g->min_weight = 0;
    g->max_weight = 1;
}

void free_graph(graph_t* g) {

    if (g->n != 0) {
        free(g->numEdges);
    }

    if (g->m != 0) {
        free(g->endV);
        free(g->edge_id);
    }

    if (g->weight_type != -1) {

        if (g->weight_type == 1) {
            free(g->int_weight_e);
        }

        if (g->weight_type == 2) {
            free(g->l_weight_e); 
        }

        if (g->weight_type == 3) {
            free(g->fl_weight_e); 
        }

        if (g->weight_type == 4) {
            free(g->dbl_weight_e); 
        }
    }

}
