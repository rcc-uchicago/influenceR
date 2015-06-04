#include "graph_partitioning.h"
#include "graph_kernels.h"
#include "graph_metrics.h"

double comm_evaluate_modularity(graph_t* g, comm_list_bc_t* comm_list, 
        attr_id_t num_components) {

    attr_id_t u, v, j;
    attr_id_t n, m;
    attr_id_t cid;
    double m2_inv; 
    attr_id_t u_edge_start, u_edge_end, u_degree;
    double mod;

    n = g->n;
    m = g->m; 
    m2_inv = 1.0/(g->m);
    mod = 0.0;

    for (j=0; j<num_components; j++) {
        comm_list[j].a = 0.0;
        comm_list[j].e = 0.0;
    }

    for(u=0; u<n; u++) {
        cid = g->cvl[u].comm_id;
        u_edge_start = g->cvl[u].num_edges;
        u_edge_end   = g->cvl[u+1].num_edges;
        u_degree     = u_edge_end - u_edge_start;

        for(j=u_edge_start; j<u_edge_end; j++) {
            v = g->cel[j].dest;
            if (g->cvl[v].comm_id == cid)
                comm_list[cid].e += 1.0;
        }
        comm_list[cid].a += (double) u_degree;
    }

    for (j=0; j<num_components; j++) {
        mod += comm_list[j].e - m2_inv * comm_list[j].a * comm_list[j].a; 
    }
    mod = mod * m2_inv;
    return mod;
}

void  remove_maxbc_edge(graph_t *g, comm_list_bc_t *comm_list, 
        attr_id_t num_components, attr_id_t num_bc_runs, 
        edge_t* ebc_edge, attr_id_t* maxbc_component) {

    attr_id_t i, j;
    double mbc_val;
    attr_id_t mbc_eid;
    attr_id_t mbc_esrc, mbc_edest;
    attr_id_t mbc_component;
    attr_id_t i_start_edge, i_end_edge;
    attr_id_t n;
    attr_id_t comm_id;
    attr_id_t prev_maxbc_component;
    n = g->n;

    /* Initialize for first run */
    if (num_bc_runs == 1) {
        for (i=0; i<n; i++) {        
            comm_id = g->cvl[i].comm_id;
            for (j=g->cvl[i].num_edges; j<g->cvl[i+1].num_edges; j++) {
                if (g->cel[j].cval > comm_list[comm_id].mbc_val) {
                    comm_list[comm_id].mbc_val = g->cel[j].cval;
                    comm_list[comm_id].mbc_esrc = i;
                    comm_list[comm_id].mbc_eid = g->cel[j].eid;
                }
            }
        }
    } else {

        prev_maxbc_component = *maxbc_component;
        /* Two components are updated in prev BC iteration, num_components-1
           and prev_maxbc_component. Update the max centrality scores of these
           components.
         */
        comm_list[prev_maxbc_component].mbc_val = -1;
        comm_list[num_components-1].mbc_val = -1;
        for (i=0; i<n; i++) {        
            comm_id = g->cvl[i].comm_id;
            if ((comm_id != prev_maxbc_component) && (comm_id !=
                        num_components-1))
                continue; 
            for (j=g->cvl[i].num_edges; j<g->cvl[i+1].num_edges; j++) {
                if (g->cel[j].cval > comm_list[comm_id].mbc_val) {
                    comm_list[comm_id].mbc_val = g->cel[j].cval;
                    comm_list[comm_id].mbc_esrc = i;
                    comm_list[comm_id].mbc_eid = g->cel[j].eid;
                }
            }
        }
    }

    /* find edge with max edge bc value */ 
    mbc_val = comm_list[0].mbc_val;
    mbc_eid = comm_list[0].mbc_eid;
    mbc_esrc = comm_list[0].mbc_esrc;
    mbc_component = 0;

    for (i=1; i<num_components; i++) {
        if (comm_list[i].mbc_val > mbc_val) {     
            mbc_val  = comm_list[i].mbc_val;
            mbc_component = i;
            mbc_eid  = comm_list[i].mbc_eid;
            mbc_esrc = comm_list[i].mbc_esrc;     
        }
    }

    assert(mbc_val != -1);

    /* Mark the corresponding edge as deleted */
    i = mbc_esrc;
    i_start_edge = g->cvl[i].num_edges;
    i_end_edge = g->cvl[i+1].num_edges;
    for (j=i_start_edge; j<i_end_edge; j++) {
        if (g->cel[j].eid == mbc_eid) {
            mbc_edest = g->cel[j].dest; 
            g->cel[j].mask = num_bc_runs;
            g->cel[j].cval = -5;
        } 
    }

    i = mbc_edest;
    i_start_edge = g->cvl[i].num_edges;
    i_end_edge = g->cvl[i+1].num_edges;
    for (j=i_start_edge; j<i_end_edge; j++) {
        if (g->cel[j].eid == mbc_eid) {
            g->cel[j].mask = num_bc_runs;
            g->cel[j].cval = -5;
        } 
    }
    /* 
       for (i=0; i<n; i++) {        
       for (j=g->cvl[i].num_edges; j<g->cvl[i+1].num_edges; j++) {
       fprintf(stderr, "%5.5lf %4d (%2d %2d) %d\n",
       g->cel[j].cval, g->cel[j].eid,
       i, g->cel[j].dest, g->cel[j].mask);
       }
       }
     */
    *maxbc_component = mbc_component;

}

/* Community detection algorithm that optimizes modularity, and is based on 
   iterative removal of edges with high betweenness in the network */
void modularity_betweenness(graph_t *g, attr_id_t *membership, 
        attr_id_t *num_communities, double *modularity, 
        double sampling_val) {

    attr_id_t n, m;
    long i, j;    
    attr_id_t num_components;
    double curr_modularity, prev_modularity;
    comm_list_bc_t* comm_list;
    edge_t ebc_edge;
    int num_bc_runs;
    attr_id_t curr_component1, curr_component2, maxbc_component;
    int split;
    attr_id_t* ebc_eval_data1;
    double* ebc_eval_data2;
    double max_modularity;
    attr_id_t max_modularity_comp_num;

    n = g->n;
    m = g->m;

    /* Initialize the graph representation we will use in this algorithm */ 
    g->cvl = (c_vert_t *) calloc(n + 1, sizeof(c_vert_t));
    g->cel = (c_edge_t *) calloc(m, sizeof(c_edge_t));
    assert(g->cvl != NULL);
    assert(g->cel != NULL);

    for (i=0; i<n; i++) {
        g->cvl[i].num_edges = g->numEdges[i];
        for (j=g->numEdges[i]; j<g->numEdges[i+1]; j++) {
            g->cel[j].dest = g->endV[j];
            g->cel[j].eid  = g->edge_id[j];
        }
    }

    g->cvl[n].num_edges = g->numEdges[n];
    /* Initially, all vertices belong to one giant community */

    /* Run connected components */
    num_components = aux_connected_components_init(g);
    fprintf(stderr, "The network has %d connected components.\n",
            num_components);

    /* We store the community splits to reconstruct the hierarchical 
       community dendrogram */
    /* Every community stores the ID of its parent, and also 
     * two variables to count the  */ 
    /* The max. number of communities is n, the number of vertices in the
     * network */
    comm_list = (comm_list_bc_t *) calloc(n, sizeof(comm_list_bc_t));
    curr_modularity = comm_evaluate_modularity(g, comm_list, num_components);
    prev_modularity = 0;
    max_modularity = 0;
    num_bc_runs = 0;

    /* Initially run betweenness computation for all connected components */
    curr_component1 = curr_component2 = -1;
    split = 0;

    /* Pre-allocate memory for the edge BC computation in every iteration */

    /* start the iterative edge-betweenness based partitioning */
    while (1) {

        evaluate_edge_centrality_bcpart(g, ebc_eval_data1, ebc_eval_data2, 
                comm_list, num_components, curr_component1, curr_component2);

        num_bc_runs++;

        remove_maxbc_edge(g, comm_list, num_components, num_bc_runs, &ebc_edge, 
                &maxbc_component);

        split = aux_connected_components_update(g, num_components,
                maxbc_component);
        if (split) {
            curr_component1 = maxbc_component;
            curr_component2 = num_components;
            comm_list[num_components].p = maxbc_component;
            num_components++;
            curr_modularity = comm_evaluate_modularity(g, comm_list, 
                    num_components);
            if (curr_modularity > max_modularity) {
                max_modularity = curr_modularity;
                max_modularity_comp_num = num_components-1;
            }
            if (curr_modularity < max_modularity - 0.25)
                break;
        }

    } 

    for (i=max_modularity_comp_num+1; i<num_components; i++) {
        j = i;
        while (j > max_modularity_comp_num) {
            j = comm_list[j].p;
        }
        comm_list[i].p = j;
    }

    for (i=0; i<n; i++) {
        if (g->cvl[i].comm_id < max_modularity_comp_num+1)
            membership[i] = g->cvl[i].comm_id;
        else
            membership[i] = comm_list[g->cvl[i].comm_id].p;
    }
    *num_communities = max_modularity_comp_num+1;
    *modularity = max_modularity;
}
