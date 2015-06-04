#if !defined(MODULARITY_GREEDY_AGGLOMERATIVE_HEADER_)
#define MODULARITY_GREEDY_AGGLOMERATIVE_HEADER_

#undef VERBOSE

#undef DEBUG

#undef DEBUG_LIJ
#undef DEBUG_SIJ
#undef DEBUG_STATS
#undef DEBUG_MB
#undef DEBUG_MAXHEAP
#undef DEBUG_MAXKEY
#undef DEBUG_LIN

#define DO_HEAP2
#define DO_ALIN
#define DO_COMM_WT

#define DO_STOPPING_LIJ

#define CULL_DEGREE_ONE
#define CULL_THRESHOLD 3

#define PREPROC_BICONN 0

#define NUMKEYTYPE 5
#define KEY_CNM   0
#define KEY_MB    1
#define KEY_RAT   2
#define KEY_MBRAT 3
#define KEY_LIN   4

#ifdef MAXDOUBLE
#define INFTY MAXDOUBLE
#else
#define INFTY 1.79769313486231570e+308
#endif

#define NO_MAX_KEY -INFTY

#if defined(__GNUC__) && __GNUC__ >= 4
#define VIS_LOCAL __attribute__ ((visibility ("hidden")))
#else
#define VIS_LOCAL
#endif

extern VIS_LOCAL int keytype; /* 0:dq This is the index to use for the value used in the heap */

typedef struct {
    int weighted;
    union {
        attr_id_t A;
        double M;
    };
    double origM;
} weight_t;

#define MAXHISTORY 20
extern VIS_LOCAL double MaxKeyHistory[MAXHISTORY];
extern VIS_LOCAL int MaxKeyHistoryIdx;
extern VIS_LOCAL int MaxKeyHistoryFilled;
extern VIS_LOCAL double sumMKH, sumMKH2;

/* Number of standard deviations to consider in the MB algorithm. */
#define NSTDDEV 1.5

/* Community adjacency structure -- sorted array */

typedef struct {
    attr_id_t comm_id; /* the list is sorted by this value */  
    union {
        attr_id_t lij;     /* Number of edges between communities i and j */
        double    Xij;     /* Sum of edge weights between communities i and j */
    };
    double value[NUMKEYTYPE];   
    /* We are using the "value" array to represent the scores related to
       each community that are potential keys for the maxheap. 
       value[0]: This is dq: the increase in modularity for this pair
       value[1]: This is dq/stddev(dq), where the stddev is with respect
       to this community's adjacent communities 
       value[2]: This is the size ratio scaled CNM dq: dq * min(s1,s2) / max(s1,s2)
       value[3]: This is the size ratio scaled dq/stddev(dq)
       value[4]: This is the linear model
       "keytype" is an index into the value[] array to use for the heap,
0:dq, 1:dq/stddev(dq) 
     */
} aggc_adjcomm_t;

typedef struct {
    attr_id_t comm_id;
    double val;
#ifdef DO_HEAP2
    double val2;
#endif
} aggc_heap_t;

/* Structure representing a community and its adjacencies */

typedef struct {
    aggc_adjcomm_t *adjcomm; /* communities adjacent to current comm */
    attr_id_t degree;         /* number of adjacent communities */
    attr_id_t max_key_idx;    /* neighbor with max key value */
    double max_key;           /* max key value */
    double a;                 /* a, sum of degrees of all vertices in community */
#ifdef DO_ALIN
    double a_lin;             /* modified a for KEY_LIN, during updates, subtracts 2*edge */
#endif
    attr_id_t max_size;       /* max size of the adjcomm array, array
                                 needs to be resized id degree = max_size */
    attr_id_t comm_size;      /* number of vertices that belong to a community */
#ifdef DO_COMM_WT
    double comm_wt;           /* Weight (sum of edges) of a community */
#endif
    attr_id_t parent_id;      /* parent of a vertex (set once community
                                 is merged) */
    attr_id_t resized;        /* flag to indicate if a community
                                 adjacency array has been resized after
                                 creation. Set to 1 on resize */
    /* char[20] filler; */    /* can align it to cache line size */
} aggc_comm_t;

/* Implicit max heap struct for storing the modularity delta values
   for community pairs */
typedef struct {
    aggc_heap_t *heap;
    attr_id_t *index;
    attr_id_t n;
#ifdef DO_HEAP2
    int use2;
#endif
} aggc_maxheap_t;

void aggc_maxheap_sift_up (aggc_maxheap_t*, attr_id_t) VIS_LOCAL;
void aggc_maxheap_check(aggc_maxheap_t*) VIS_LOCAL;
double aggc_merge_communities(aggc_comm_t*, aggc_maxheap_t*, aggc_adjcomm_t*,
			      attr_id_t, weight_t *) VIS_LOCAL;
attr_id_t aggc_compute_membership(aggc_comm_t*, attr_id_t*, attr_id_t) VIS_LOCAL;

/* Shared statically to encourage compilers to be smarter... */
static int
comp_comm(const void *v1, const void *v2)
{
    /* This function is used to by qsort() to sort the array of
       aggc_adjcomm_t structures */

    const aggc_adjcomm_t *p1, *p2;
    p1 = (aggc_adjcomm_t *)v1;
    p2 = (aggc_adjcomm_t *)v2;

    if (p1->comm_id == p2->comm_id)
	return 0;
    else
	return ((p1->comm_id < p2->comm_id) ? -1 : 1);
}

#endif /* MODULARITY_GREEDY_AGGLOMERATIVE_HEADER_ */
