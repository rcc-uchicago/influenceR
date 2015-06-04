#ifndef _GRAPH_DEFS_H
#define _GRAPH_DEFS_H

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#define OMP(x) _Pragma(x)
#else
#define OMP(x)
#endif

#if ENABLE_64BIT_VID
    typedef long attr_id_t;
#else
    typedef int attr_id_t;
#endif

#define ARRAY_INIT_SIZE 8

#if !HAVE_LOG2
#ifndef log2
#define log2(d) (log(d)/log(2.0))
#endif
#endif

typedef struct {
    attr_id_t* vals;
    attr_id_t* aux_info;
    long count;
    long max_size;
} dyn_array_t;

typedef struct {
    /* value is the key */
    attr_id_t *vals;
    attr_id_t count;
    attr_id_t max_size;
} adj_bheap_t;

void adj_bheap_insert(adj_bheap_t* h, attr_id_t val);


void dyn_array_insert(dyn_array_t* A, attr_id_t val);
void dyn_array_delete(dyn_array_t* A, attr_id_t val);
void dyn_array_init(dyn_array_t* A);
void dyn_array_clear(dyn_array_t* A);
void dyn_array_free(dyn_array_t* A);

void sorted_dyn_array_insert(dyn_array_t* A, attr_id_t val);
void sorted_dyn_array_delete(dyn_array_t* A, attr_id_t val);
void sorted_dyn_array_init(dyn_array_t* A);
void sorted_dyn_array_clear(dyn_array_t* A);
void sorted_dyn_array_free(dyn_array_t* A);

typedef struct {
    attr_id_t dest;
    attr_id_t eid;
    attr_id_t mask;
    double cval;
} c_edge_t;

typedef struct {
    attr_id_t num_edges;
    attr_id_t comm_id;
} c_vert_t;

typedef struct {
    attr_id_t p;
    double a;
    double e;
    attr_id_t mbc_eid;
    double mbc_val;
    attr_id_t mbc_esrc;
} comm_list_bc_t;


typedef struct
{
    attr_id_t src;
    attr_id_t dest;
    attr_id_t eid;
} edge_t;

typedef struct
{
    dyn_array_t* neighbors;
    attr_id_t degree;
} vert_t;

typedef struct {
    attr_id_t* list;
    attr_id_t count;
    attr_id_t degree;
} plist_t;

typedef struct struct_node
{
    int id;
    int aux_val;
    struct struct_node *next;
} node_t;

typedef struct
{
    node_t *head;
    node_t *tail;
    int size;
} list_t;


/* The graph data structure*/
typedef struct
{
    /***
     The minimal graph repesentation consists of:
     n        -- the number of vertices
     m        -- the number of edges
     endV     -- an array of size 2*m (m for directed graphs) that stores the 
                 destination ID of an edge <src->dest>.
     numEdges -- an array of size n+1 that stores the degree 
                 (out-degree in case of directed graphs) and pointers to
                 the endV array. The degree of vertex i is given by 
                 numEdges[i+1]-numEdges[i], and the edges out of i are
                 stored in the contiguous block endV[numEdges[i] ..
                 numEdges[i+1]].
     Vertices are ordered from 0 in our internal representation
     ***/
    long n;
    long m;
    attr_id_t* endV;
    attr_id_t* numEdges;
    
    int undirected;
    int zero_indexed;

    /***
     Other useful constructs that can be initialized when needed
     ***/
    edge_t* elist;
    edge_t* elist_aux;
    vert_t* vlist;
    vert_t* vlist_aux;

    /* An internal id for each edge, useful in case of undirected networks */
    attr_id_t* edge_id;

    /* For directed graphs, endBackV can be used to store edges pointing into a
       vertex, and numBackEdges the in-degree */
    attr_id_t* endBackV;
    attr_id_t* numBackEdges;
   
    /* Data representation used in some centrality and community identification
     * routines */
    c_vert_t* cvl;
    c_edge_t* cel;

    /* Vertex weights */
    int* int_weight_v;
    float* fl_weight_v;
    double* dbl_weight_v;
    long* l_weight_v;

    /* Edge weights */
    int weight_type;
    int* int_weight_e;
    float* fl_weight_e;
    double* dbl_weight_e;
    long* l_weight_e;

    double min_weight;
    double max_weight;

    /* Fine-grained locking support, currently using 
       OpenMP mutex locks */
#ifdef _OPENMP
    omp_lock_t* vLock;
    omp_lock_t* eLock;
#endif

} graph_t;

/* A simple data structure to experiment with
   dynamic networks
 */
typedef struct
{
    /***
     The minimal dynamic graph repesentation consists of:
     n        -- the number of vertices
     m        -- the number of edges
     adj      -- representation of adjacencies of every vertex 
     spanf    -- a forest of rooted spanning trees for answering connectivity
     queries in the dynamic network
     num_trees-- the number of spanning trees in the forest
     Vertices are ordered from 0 in our internal representation
    ***/
    long n;
    long m;
    dyn_array_t* adj;
    adj_bheap_t* Hadj;
    attr_id_t* degree_thresh;
    attr_id_t *pred;
    attr_id_t *ts;
} dyn_graph_t;

void par_gen_RMAT_edges(dyn_graph_t* G, double* params, attr_id_t *src, attr_id_t* dest,
        attr_id_t* degree);
void dyn_ds_init_base(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree);
void dyn_ds_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree);
void dyn_ds_init_malloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest);
void dyn_vpart_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree);
void dyn_epart_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree);
void dyn_heap_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const
        attr_id_t* dest, const attr_id_t* degree);
void dyn_hybrid_init_nomalloc(dyn_graph_t* G, const attr_id_t* src, const
        attr_id_t* dest, const attr_id_t* degree);
void dyn_gr_init(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree);
void dyn_hybrid_gr_init(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree);
void dyn_gr_del(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree, attr_id_t num_dels);
void dyn_hybrid_gr_del(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        const attr_id_t* degree, attr_id_t numdels);
void dyn_st_forest_construct(dyn_graph_t* G);
void dyn_st_queries(dyn_graph_t* G, const attr_id_t* src, const attr_id_t* dest,
        attr_id_t num_dels);
void dyn_induced_subgraphs(dyn_graph_t* G, attr_id_t ts1, attr_id_t
        ts2);
void dyn_gr_traversal(dyn_graph_t* G, attr_id_t ts1, attr_id_t ts2);

#endif
