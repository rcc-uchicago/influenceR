/*
 AUTHOR: D.A. Bader and K. Madduri
 LICENSE: GPLv2 
 See: http://snap-graph.sourceforge.net/
  
 SNAP graph data structure, taken from snap-graph project.
*/

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



#endif
