#include "graph_defs.h"
#include "graph_kernels.h"

int aux_connected_components_init(graph_t* g) {

    attr_id_t num_components;
    attr_id_t i, j, u, v;
    attr_id_t n;
    attr_id_t *S;
    attr_id_t vis_current, vis_count;
    attr_id_t i_start_edge, i_end_edge;

    num_components = 0;
    n = g->n;

    S = (attr_id_t *) malloc(n * sizeof(attr_id_t));
    assert(S != NULL);

    for (i=0; i<n; i++) {
        g->cvl[i].comm_id = -1;
    }

    for (i=0; i<n; i++) {
        if (g->cvl[i].comm_id != -1) 
            continue;
         
        /* Set component ID of all vertices reachable from i */
        S[0] = i;
        vis_count = 1;
        vis_current = 0;
        g->cvl[i].comm_id = num_components;
        while (vis_current != vis_count) {
            u = S[vis_current];
            i_start_edge = g->cvl[u].num_edges;
            i_end_edge   = g->cvl[u+1].num_edges;
            for (j=i_start_edge; j<i_end_edge; j++) {
                v = g->cel[j].dest;
                if (g->cvl[v].comm_id == -1) {
                    g->cvl[v].comm_id = num_components;
                    S[vis_count++] = v;
                }
            }
            vis_current++;
        }
        num_components++;
    }
    free(S);
    return num_components;
}

int aux_connected_components_update(graph_t* g, attr_id_t num_components, 
        attr_id_t maxbc_component) {

    int split;
    attr_id_t i, j, v, u;
    attr_id_t n;
    attr_id_t *S;
    attr_id_t vis_current, vis_count;
    attr_id_t i_start_edge, i_end_edge;

    split = -1;
    n = g->n;

    S = (attr_id_t *) malloc(n * sizeof(attr_id_t));
    assert(S != NULL);

    for (i=0; i<n; i++) {
        if (g->cvl[i].comm_id != maxbc_component)
            continue;
        
        if (split == -1)  {
            split = 0;
            /* Set component ID of all vertices reachable from i */
            S[0] = i;
            vis_count = 1;
            vis_current = 0;
            g->cvl[i].comm_id = num_components;
            while (vis_current != vis_count) {
                u = S[vis_current];
                i_start_edge = g->cvl[u].num_edges;
                i_end_edge   = g->cvl[u+1].num_edges;
                for (j=i_start_edge; j<i_end_edge; j++) {
                    if (g->cel[j].mask == 0) {
                        v = g->cel[j].dest;
                        if (g->cvl[v].comm_id == maxbc_component) {
                            g->cvl[v].comm_id = num_components;
                            S[vis_count++] = v;
                        }
                    }
                }
                vis_current++;
            }
        } else if (split == 0) {
            split = 1;
            break;
        }
    }

    /* If there is no new component created, reset the component 
     * ID of vertices in S */
    if (split == 0) {
        for (i=0; i<vis_count; i++) {
            g->cvl[S[i]].comm_id = maxbc_component;
        }
    }

    free(S);
    return split;
}


void connected_components(graph_t* g, attr_id_t* component_num) {

}


