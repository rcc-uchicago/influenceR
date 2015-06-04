#include <string.h>
#include "graph_partitioning.h"
#include "graph_kernels.h"

#include "modularity_greedy_agglomerative.h"

#if defined(REALLY_BROKEN_ARITHMETIC_THAT_NEEDS_FIXED)
#define MACHINE_EPS 0.0000001
#define fpEQ(a,b)             (fabs((a)-(b))<MACHINE_EPS)
#define fpLT(a,b) (((a)<(b))&&(fabs((a)-(b))>MACHINE_EPS))
#define fpGT(a,b) (((a)>(b))&&(fabs((a)-(b))>MACHINE_EPS))
#else
#define fpEQ(a,b) ((a) == (b))
#define fpLT(a,b) ((a)<(b))
#define fpGT(a,b) ((a)>(b))
#endif

static void
aggc_print_communities(aggc_comm_t *communities, attr_id_t n, weight_t L)
{
    attr_id_t i, j;

    fprintf(stderr, "\n");
    for (i=0; i<n; i++) {
	if (communities[i].max_size > 0) {
	    fprintf(stderr, "Comm %9d, size %9d  max size %9d, degree %9d\n", i,
		    communities[i].comm_size, communities[i].max_size, communities[i].degree);
	    fprintf(stderr, "max_key %9.6f, idx %d, parent id %d, a %9.6f, resized %d\n",
		    (communities[i].max_key == NO_MAX_KEY ? -999 : communities[i].max_key),
		    communities[i].max_key_idx, communities[i].parent_id,
		    communities[i].a, communities[i].resized);
	    if (communities[i].degree > 0)
		fprintf(stderr,"[%9s %15s %15s %15s %15s %15s]\n","comm_id","Xij","value[CNM]","L*val[CNM]","val[1]","val[KEY");
	    for (j=0; j<communities[i].degree; j++) {
		fprintf(stderr, "[%9d ", communities[i].adjcomm[j].comm_id);
		if (L.weighted)
		    fprintf(stderr, "%15.9f ",  communities[i].adjcomm[j].Xij);
		else
		    fprintf(stderr, "%9d ",  communities[i].adjcomm[j].lij);
		fprintf(stderr, "%15.9f ", communities[i].adjcomm[j].value[KEY_CNM]);
#if 1
		fprintf(stderr, "%15.9f ", (L.weighted ? L.M : (double)L.A)*communities[i].adjcomm[j].value[KEY_CNM]);
		fprintf(stderr, "%15.9f ", communities[i].adjcomm[j].value[1]);
#endif
		fprintf(stderr, "%15.9f ", communities[i].adjcomm[j].value[keytype]);
		fprintf(stderr,"]\n");
	    }
	    fprintf(stderr, "\n\n");
	}
    }
    return;
}

static void
aggc_print_final_communities(aggc_comm_t *communities, aggc_maxheap_t *maxheap, weight_t L)
{
    attr_id_t i, j, k;

    fprintf(stderr, "\n");


    for (i=0; i<maxheap->n; i++) {
	k = maxheap->heap[i].comm_id;

	if (communities[k].max_size > 0) {
	    fprintf(stderr, "Comm %9d, size %9d  degree %9d\n",
		    k, communities[k].comm_size, communities[k].degree);
	    for (j=0; j<communities[k].degree; j++) {
		fprintf(stderr, "[%9d ", communities[k].adjcomm[j].comm_id);
		if (L.weighted)
		    fprintf(stderr, "%15.9f ",  communities[k].adjcomm[j].Xij);
		else
		    fprintf(stderr, "%9d ",  communities[k].adjcomm[j].lij);
		fprintf(stderr, "%15.9f ", communities[k].adjcomm[j].value[KEY_CNM]);
#if 1
		fprintf(stderr, "%15.9f ", (L.weighted ? L.M : (double)L.A)*communities[k].adjcomm[j].value[KEY_CNM]);
		fprintf(stderr, "%15.9f ", communities[k].adjcomm[j].value[1]);
#endif
		fprintf(stderr,"]\n");
	    }
	    fprintf(stderr, "\n\n");
	}
    }

    return;
}

void
aggc_maxheap_check(aggc_maxheap_t *maxheap)
{
    /* This function checks that the heap property is not violated. */
    attr_id_t i, n;
    aggc_heap_t *heap;
    int broken=0;

#ifdef DEBUG
    fprintf(stderr,"Checking heap.\n");
#ifdef DO_HEAP2
    if (maxheap->use2)
	fprintf(stderr,"Checking heap with use2.\n");
#endif
#endif

    heap    = maxheap->heap;
    n       = maxheap->n;

    for (i=0; i<n/2; i++) {
	if (2*i+1<n) {
	    if ( (fpLT(heap[i].val,heap[2*i+1].val))
#ifdef DO_HEAP2
		    || ((maxheap->use2)&&(fpEQ(heap[i].val,heap[2*i+1].val))&&(fpLT(heap[i].val2,heap[2*i+1].val2)))
#endif
	       ) {
		fprintf(stderr, "ERROR: heap property violated (n: %d i: %d  2i+1: %d) %9.6f < %9.6f\n",
			n, i, 2*i+1, heap[i].val,heap[2*i+1].val);
#ifdef DO_HEAP2
		if (maxheap->use2)
		    fprintf(stderr, "ERROR: heap property violated  val2 %9.6f < %9.6f\n",
			    heap[i].val2,heap[2*i+1].val2);
#endif
		broken=1;
	    }
	}
	if (2*i+2<n) {
	    if ( (fpLT(heap[i].val,heap[2*i+2].val))
#ifdef DO_HEAP2
		    || ((maxheap->use2)&&(fpEQ(heap[i].val,heap[2*i+2].val))&&(fpLT(heap[i].val2,heap[2*i+2].val2)))
#endif
	       ) {
		fprintf(stderr, "ERROR: heap property violated (n: %d i: %d  2i+2: %d) %9.6f < %9.6f\n",
			n, i, 2*i+2, heap[i].val,heap[2*i+2].val);
#ifdef DO_HEAP2
		if (maxheap->use2)
		    fprintf(stderr, "ERROR: heap property violated  val2 %9.6f < %9.6f\n",
			    heap[i].val2,heap[2*i+2].val2);
#endif
		broken=1;
	    }
	}
    }

    if (broken) {
	fprintf(stderr,"ERROR: Heap broken. Dumping heap.\n");
	for (i=0; i<n; i++) {
	    fprintf(stderr, "heap[%4d].val (%9.6f) "
#ifdef DO_HEAP2
		    "val2 (%9.6f) "
#endif
		    "comm_id (%4d)\n",
		    i, heap[i].val,
#ifdef DO_HEAP2
		    heap[i].val2,
#endif
		    heap[i].comm_id);
	}

	exit(-1);
    }

#ifdef DEBUG
    for (i=0 ; i<n ; i++)
	fprintf(stderr,"DEBUG maxheap val[%5d]: %9.6f\n",i,heap[i].val);
#endif

    return;
}

static void
aggc_maxkey_comm_check(aggc_comm_t *communities, attr_id_t comm1)
{

    attr_id_t i;

    for (i=0 ; i<communities[comm1].degree ; i++)
	if (communities[comm1].adjcomm[i].value[keytype] > communities[comm1].max_key)
	    fprintf(stderr,"ERROR: maxkey: %f  larger: %f (%2d, %2d)\n",
		    communities[comm1].max_key,
		    communities[comm1].adjcomm[i].value[keytype], comm1, i);
    return;
}

static void
aggc_maxkey_check(aggc_comm_t *communities, aggc_maxheap_t *maxheap)
{
    attr_id_t i, j, k;
    double x, max_val;

    max_val = maxheap->heap[0].val;

    for (i=0; i<maxheap->n; i++) {
	k = maxheap->heap[i].comm_id;
	for (j=0 ; j<communities[k].degree ; j++) {
	    x = communities[k].adjcomm[j].value[keytype];
	    if (x > max_val)
		fprintf(stderr,"ERROR: maxkey_check violated (%7.4f at (%2d, %2d, %2d)), heap[0] %7.4f\n",
			x, i, j, k, max_val);
	}
    }

    return;
}

void
aggc_maxheap_sift_up(aggc_maxheap_t *maxheap, attr_id_t idx)
{
    /* When a heap key value increases, this function will start at that
       key at position "idx" and sift up by performing swaps moving
       up towards the top node of the heap. This ensures that the heap
       property is maintained. */

    aggc_heap_t *heap;
    attr_id_t *index;
    attr_id_t root, parent, u, v, n, pos;
    double tmp;
#ifdef DO_HEAP2
    double tmp2;
#endif

    heap    = maxheap->heap;
    index   = maxheap->index;
    n       = maxheap->n;

    root = idx;

#ifdef DEBUG
    fprintf(stderr,"DEBUG SIFTUP: idx: %d  vals (%f, %f)\n",idx,heap[idx].val, heap[idx].val2);
#endif

#ifdef DO_HEAP2
    if (!maxheap->use2) {
#endif
	while (root > 0) {
	    parent = (root-1)/2;
	    if (heap[parent].val < heap[root].val) {
		/* swap */
		tmp = heap[root].val;
		v   = heap[root].comm_id;
		u   = heap[parent].comm_id;

		heap[root].val = heap[parent].val;
		heap[root].comm_id = u;

		heap[parent].val = tmp;
		heap[parent].comm_id = v;

		pos = index[v];
		index[v] = index[u];
		index[u] = pos;

		root = parent;

	    } else break;
	}
#ifdef DO_HEAP2
    }
    else {
	while (root > 0) {
	    parent = (root-1)/2;
	    if ( (fpLT(heap[parent].val,heap[root].val)) ||
		    ((fpEQ(heap[parent].val,heap[root].val)) && (fpLT(heap[parent].val2,heap[root].val2)))
	       )  {
		/* swap */
		tmp  = heap[root].val;
		tmp2 = heap[root].val2;
		v   = heap[root].comm_id;
		u   = heap[parent].comm_id;

		heap[root].val  = heap[parent].val;
		heap[root].val2 = heap[parent].val2;
		heap[root].comm_id = u;

		heap[parent].val  = tmp;
		heap[parent].val2 = tmp2;
		heap[parent].comm_id = v;

		pos = index[v];
		index[v] = index[u];
		index[u] = pos;

		root = parent;

	    } else break;
	}
    }
#endif
    return;
}

void
aggc_maxheap_sift_down(aggc_maxheap_t *maxheap, attr_id_t idx)
{
    /* When a heap key value decreases, this function will start at that
       key at position "idx" and sift down by performing swaps moving
       down towards the leaves of the heap. This ensures that the heap
       property is maintained. */

    aggc_heap_t *heap;
    attr_id_t *index;
    attr_id_t root, child, u, v, n, pos;
    double tmp;
#ifdef DO_HEAP2
    double tmp2;
#endif

    heap    = maxheap->heap;
    index   = maxheap->index;
    n       = maxheap->n;

    root = idx;

#ifdef DEBUG
    fprintf(stderr,"DEBUG SIFTDN: idx: %d  vals (%f, %f)\n",idx,heap[idx].val, heap[idx].val2);
#endif

#ifdef DO_HEAP2
    if (!maxheap->use2) {
#endif
	while (2*root+1 < n) {
	    child = 2*root+1;
	    if (child+1 < n)
		if (heap[child].val < heap[child+1].val)
		    child++;
	    if (heap[root].val < heap[child].val) {
		/* swap */
		tmp = heap[root].val;
		v  = heap[root].comm_id;
		u  = heap[child].comm_id;

		heap[root].val = heap[child].val;
		heap[root].comm_id = u;

		heap[child].val = tmp;
		heap[child].comm_id = v;

		pos = index[v];
		index[v] = index[u];
		index[u] = pos;

		root = child;

	    } else break;
	}
#ifdef DO_HEAP2
    }
    else {
	while (2*root+1 < n) {
	    child = 2*root+1;
	    if (child+1 < n)
		if ( (fpLT(heap[child].val,heap[child+1].val)) ||
			((fpEQ(heap[child].val, heap[child+1].val)) && (fpLT(heap[child].val2,heap[child+1].val2)))
		   )
		    child++;
	    if ((fpLT(heap[root].val,heap[child].val)) ||
		    ((fpEQ(heap[root].val, heap[child].val)) && (fpLT(heap[root].val2,heap[child].val2)))
	       ) {
		/* swap */
		tmp  = heap[root].val;
		tmp2 = heap[root].val2;
		v  = heap[root].comm_id;
		u  = heap[child].comm_id;

		heap[root].val  = heap[child].val;
		heap[root].val2 = heap[child].val2;
		heap[root].comm_id = u;

		heap[child].val  = tmp;
		heap[child].val2 = tmp2;
		heap[child].comm_id = v;

		pos = index[v];
		index[v] = index[u];
		index[u] = pos;

		root = child;

	    } else break;
	}
    }
#endif
    return;
}

static void
aggc_maxheap_remove(aggc_maxheap_t *maxheap, attr_id_t comm_id)
{
    /* When a community with id comm_id is merged, this function removes
       it from the heap. */
    aggc_heap_t *heap;
    attr_id_t *index;
    attr_id_t idx, new_comm_id, n;
    double old_val;
#ifdef DO_HEAP2
    double old_val2;
#endif

    heap    = maxheap->heap;
    index   = maxheap->index;
    n       = maxheap->n;

    idx = index[comm_id];
#if 0
    assert(idx != -1);
#else
    if (idx < 0) return;
#endif

#ifdef DO_HEAP2
    if (!maxheap->use2) {
#endif
	old_val = heap[idx].val;

	/* Move last value to this position */
	new_comm_id = heap[n-1].comm_id;
	heap[idx].comm_id = new_comm_id;
	heap[idx].val     = heap[n-1].val;

	index[new_comm_id] = idx;
	index[comm_id]     = -1;

	maxheap->n = n-1;

	if (fpGT(heap[idx].val,old_val))
	    aggc_maxheap_sift_up(maxheap, idx);
	else
	    aggc_maxheap_sift_down(maxheap, idx);

#ifdef DEBUG_MAXHEAP
	fprintf(stderr,"DEBUG_MAXHEAP Checking heap: maxheap_remove end\n");
	aggc_maxheap_check(maxheap);
#endif

#ifdef DO_HEAP2
    }
    else {
	old_val  = heap[idx].val;
	old_val2 = heap[idx].val2;

	/* Move last value to this position */
	new_comm_id = heap[n-1].comm_id;
	heap[idx].comm_id = new_comm_id;
	heap[idx].val     = heap[n-1].val;
	heap[idx].val2    = heap[n-1].val2;

	index[new_comm_id] = idx;
	index[comm_id]     = -1;

	maxheap->n = n-1;

	if ((fpGT(heap[idx].val,old_val)) ||
		((fpEQ(heap[idx].val, old_val)) && (fpGT(heap[idx].val2,old_val2))) )
	    aggc_maxheap_sift_up(maxheap, idx);
	else
	    aggc_maxheap_sift_down(maxheap, idx);
    }
#endif

    return;
}

static attr_id_t
aggc_adjcomm_pos(aggc_adjcomm_t *adjcomm, attr_id_t val, attr_id_t degree)
{
    /* Given an attr_id_t of val, use binary search to find its position
       in the sorted adjacent community list adjcomm. */
    attr_id_t low, high, mid;

#if 0
    fprintf(stderr,"POS Looking for: %2d  in list of size %2d\n",val,degree);
#endif

    low = 0;
    high = degree-1;
    while (low < high) {
	mid = low + ((high-low)/2);
	if (adjcomm[mid].comm_id < val)
	    low = mid + 1;
	else
	    high = mid;
    }

#if 0
    {
	attr_id_t i;
	high = -INFTY;
	for (i=0 ; i<degree ; i++) {
	    fprintf(stderr,"POS  %2d %2d\n",i,adjcomm[i].comm_id);
	    if (adjcomm[i].comm_id == val) high=i;
	}
	assert (high==low);
    }
#endif

    return low;
}

static void
aggc_adjcomm_resize(aggc_comm_t *communities, attr_id_t comm_id, attr_id_t incr)
{
    /* When we merge a smaller community into a larger one, the list of
       adjacent communities may need a larger array. This is a resize
       procedure that increases the size by the amount "incr". */
    attr_id_t new_size;
    aggc_adjcomm_t *new_adjcomm, *old_adjcomm;

    old_adjcomm = communities[comm_id].adjcomm;

    new_size = communities[comm_id].max_size + 2*incr;

    communities[comm_id].max_size = new_size;

    new_adjcomm = (aggc_adjcomm_t *)malloc(new_size * sizeof(aggc_adjcomm_t));
    assert(new_adjcomm != NULL);
    memcpy(new_adjcomm, old_adjcomm, communities[comm_id].degree * sizeof(aggc_adjcomm_t));
    if (communities[comm_id].resized == 1)
	free(old_adjcomm);
    else
	communities[comm_id].resized = 1;
    communities[comm_id].adjcomm = new_adjcomm;

    return;
}

static void
aggc_update_values(aggc_comm_t *communities, attr_id_t comm1, attr_id_t idx,
		   attr_id_t comm2, weight_t L, attr_id_t n)
{
    double lhat, rij, cratio;
    aggc_adjcomm_t *comm1_adjcomm;
    attr_id_t comm1_size, comm2_size;

    /* idx is the comm1 index for comm2 */

    if (keytype==KEY_CNM) return;

    comm1_adjcomm = communities[comm1].adjcomm;
    comm1_size = communities[comm1].comm_size;
#if 1
    assert(comm2 == comm1_adjcomm[idx].comm_id);
#endif
    comm2_size = communities[comm2].comm_size;

#ifdef DEBUG_MB
    fprintf(stderr,"DEBUG_MB update %2d %2d values comm1_size: %2d  comm2_size: %2d\n",
	    comm1, comm2, comm1_size, comm2_size);
#endif

    switch (keytype) {
	case KEY_MB:
	case KEY_RAT:
	case KEY_MBRAT:

	    /* value[1] update */

	    /* This is variance */
	    rij = 0.5 * (L.weighted ? L.M : (double)L.A) * comm1_adjcomm[idx].value[KEY_CNM];
	    lhat = (L.weighted ? comm1_adjcomm[idx].Xij : (double)comm1_adjcomm[idx].lij) - rij;

#ifdef DEBUG_LIJ
	    if (L.weighted)
		fprintf(stderr,"DEBUG_LIJ update %2d %2d values lij: %f L:%2d qij: %9.6f .. Rij: %7.4f  lij-Rij: %7.4f sqrt(lijH): %7.4f\n",
			comm1, comm1_adjcomm[idx].comm_id, comm1_adjcomm[idx].lij, L, comm1_adjcomm[idx].value[KEY_CNM],
			rij, lhat, sqrt(lhat));
	    else
		fprintf(stderr,"DEBUG_LIJ update %2d %2d values lij: %2d L:%2d qij: %9.6f .. Rij: %7.4f  lij-Rij: %7.4f sqrt(lijH): %7.4f\n",
			comm1, comm1_adjcomm[idx].comm_id, comm1_adjcomm[idx].lij, L, comm1_adjcomm[idx].value[KEY_CNM],
			rij, lhat, sqrt(lhat));
#endif

	    if (lhat<0) {
		fprintf(stderr,"Error: variance=lij-Rij (%f) < 0, exiting.\n", lhat);
		exit(-1);
	    }

	    assert(lhat >= 0);

	    lhat = sqrt(lhat);

	    if (lhat < 0.000001) fprintf(stderr,"Warning: update values lhat < 0.000001\n");

	    comm1_adjcomm[idx].value[KEY_MB] =  rij / lhat;

	    /* value[KEY_RAT] update */

	    if (comm1_size < comm2_size)
		cratio = (double)comm1_size/(double)comm2_size;
	    else
		cratio = (double)comm2_size/(double)comm1_size;
	    comm1_adjcomm[idx].value[KEY_RAT]   =  comm1_adjcomm[idx].value[KEY_CNM] * cratio;
	    /* value[KEY_MBRAT] update */
	    comm1_adjcomm[idx].value[KEY_MBRAT] =  comm1_adjcomm[idx].value[KEY_MB]  * cratio;

	    break;
	case KEY_LIN:
#ifdef DO_ALIN
	    comm1_adjcomm[idx].value[KEY_LIN] =
		2.0 *
		(comm1_adjcomm[idx].Xij - (L.origM * (communities[comm1].a_lin + communities[comm2].a_lin)/(double)n) + (L.M /(double)(n*n)) );
#else
	    comm1_adjcomm[idx].value[KEY_LIN] =
		2.0 *
		(comm1_adjcomm[idx].Xij - (L.origM * (communities[comm1].a     + communities[comm2].a    )/(double)n) + (L.M /(double)(n*n)) );
#endif

#ifdef DEBUG_LIN
	    fprintf(stderr,"DEBUG_LIN Init up %4d %4d  L (%f) Xij: %15f i.a: %8f j.a: %8f n: %d n2: %d v: %15f \n",
		    comm1, comm2,
		    L.M,
		    comm1_adjcomm[idx].Xij,
		    communities[comm1].a_lin,
		    communities[comm2].a_lin,
		    n,
		    n*n,
		    comm1_adjcomm[idx].value[KEY_LIN]);
#endif
	    break;
	case KEY_CNM:
	default:
	    fprintf(stderr,"ERROR; keytype unknown\n");

    }


    return;
}

static void
aggc_update_dq_p1(aggc_comm_t *communities, aggc_maxheap_t *maxheap,
		  attr_id_t comm1, attr_id_t comm2, weight_t L,
		  double new_dq_val, weight_t new_lij)
{
    /* Update dq value of adjcomm (comm1 -> comm2) with value new_dq_val */

    attr_id_t i, comm1_max_key_idx, comm1_max_key_adjcomm_id, pos, comm1_degree, *index, heap_pos;
    aggc_adjcomm_t *comm1_adjcomm;
    aggc_heap_t *heap;
    double comm1_max_key, comm1_max_key_old, key_val;
#ifdef DO_HEAP2
    double comm1_max_key_old2;
#endif

#ifdef DEBUG_MB
    fprintf(stderr,"DEBUG_MB UpdateP1 %2d %2d\n",comm1, comm2);
#endif

    heap  = maxheap->heap;
    index = maxheap->index;

    heap_pos = index[comm1];
    if (heap_pos < 0) return;

    comm1_max_key_idx = communities[comm1].max_key_idx;
    comm1_max_key     = communities[comm1].max_key;
    comm1_max_key_old = comm1_max_key;
#ifdef DO_HEAP2
    comm1_max_key_old2 = communities[comm1].adjcomm[comm1_max_key_idx].value[0];
#endif

    assert(comm1_max_key_idx >= 0);
    comm1_adjcomm = communities[comm1].adjcomm;
    comm1_degree  = communities[comm1].degree;
    comm1_max_key_adjcomm_id = comm1_adjcomm[comm1_max_key_idx].comm_id;

    if (comm2 == comm1_max_key_adjcomm_id) {
	comm1_adjcomm[comm1_max_key_idx].value[KEY_CNM] = new_dq_val;
	if (new_lij.weighted)
	    comm1_adjcomm[comm1_max_key_idx].Xij      = new_lij.M;
	else
	    comm1_adjcomm[comm1_max_key_idx].lij      = new_lij.A;

#ifdef DEBUG_LIJ
	fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (update p1 comm2=comm1_max)\n",comm1,comm1_adjcomm[comm1_max_key_idx].comm_id, new_lij);
#endif

#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB bef UV %2d %2d p1 c1_max_key_idx: %2d\n",
		comm1, comm2, comm1_max_key_idx);
#endif
	aggc_update_values(communities, comm1, comm1_max_key_idx, comm2, L, maxheap->n - 1);

	if (comm1_adjcomm[comm1_max_key_idx].value[keytype] > comm1_max_key) {
	    comm1_max_key = comm1_adjcomm[comm1_max_key_idx].value[keytype];
	    communities[comm1].max_key = comm1_adjcomm[comm1_max_key_idx].value[keytype];
	}
	else {
	    comm1_max_key     = NO_MAX_KEY;
	    comm1_max_key_idx = -1;
	    /* Rescan adjacencies to determine largest key value */
	    for (i=0; i<comm1_degree; i++) {
		key_val = comm1_adjcomm[i].value[keytype];
		if (key_val > comm1_max_key) {
		    comm1_max_key_idx = i;
		    comm1_max_key     = key_val;
		}
	    }
	    communities[comm1].max_key     = comm1_max_key;
	    communities[comm1].max_key_idx = comm1_max_key_idx;
	}

	/* Sift up/down value in max heap */
	heap_pos = index[comm1];
	heap[heap_pos].val = comm1_max_key;

#ifdef DO_HEAP2
	if (maxheap->use2)
	    heap[heap_pos].val2 = comm1_adjcomm[comm1_max_key_idx].value[0];
#endif

	if ( (fpGT(comm1_max_key,comm1_max_key_old))
#ifdef DO_HEAP2
		|| ((maxheap->use2)&&(fpEQ(comm1_max_key,comm1_max_key_old))&&(fpGT(heap[heap_pos].val2,comm1_max_key_old2)))
#endif
	   )
	    aggc_maxheap_sift_up(maxheap, heap_pos);
	else
	    aggc_maxheap_sift_down(maxheap, heap_pos);

#ifdef DEBUG_MAXHEAP
	fprintf(stderr,"DEBUG_MAXHEAP Checking heap: p1 A\n");
	aggc_maxheap_check(maxheap);
#endif

    }
    else {
	/* Locate this adjacency */
	pos = aggc_adjcomm_pos(comm1_adjcomm, comm2, comm1_degree);

	/* Update its key value */
	comm1_adjcomm[pos].value[KEY_CNM] = new_dq_val;
	if (new_lij.weighted)
	    comm1_adjcomm[pos].Xij      = new_lij.M;
	else
	    comm1_adjcomm[pos].lij      = new_lij.A;

#ifdef DEBUG_LIJ
	fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (update p1-else)\n",comm1,comm1_adjcomm[pos].comm_id, new_lij);
#endif

#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB bef UV %2d %2d p1-else pos: %2d\n",
		comm1, comm2, pos);
#endif
	aggc_update_values(communities, comm1, pos, comm2, L, maxheap->n - 1);

	if (comm1_adjcomm[pos].value[keytype] > comm1_max_key) {
	    /* Update global adj value */
	    communities[comm1].max_key = comm1_adjcomm[pos].value[keytype];
	    communities[comm1].max_key_idx = pos;

	    /* Sift up/down value in max heap */
	    heap[heap_pos].val = communities[comm1].max_key;
#ifdef DO_HEAP2
	    if (maxheap->use2)
		heap[heap_pos].val2 = comm1_adjcomm[pos].value[0];
#endif
	    aggc_maxheap_sift_up(maxheap, heap_pos);

#ifdef DEBUG_MAXHEAP
	    fprintf(stderr,"DEBUG_MAXHEAP Checking heap: p1 B\n");
	    aggc_maxheap_check(maxheap);
#endif

	}
    }

    return;
}


static void
aggc_update_dq_p2(aggc_comm_t *communities, aggc_maxheap_t* maxheap,
		  attr_id_t comm1, attr_id_t comm2, attr_id_t new_comm_id, weight_t L,
		  double new_dq_val, weight_t new_lij)
{
    /* Update dq value of adjcomm (comm1 -> comm2 -> new_comm_id) with value new_dq_val */

    attr_id_t i, comm1_max_key_idx, comm1_max_key_adjcomm_id, pos, comm1_degree, *index, heap_pos;
    aggc_adjcomm_t *comm1_adjcomm, new_adjcomm;
    aggc_heap_t *heap;
    double comm1_max_key, comm1_max_key_old, key_val;
#ifdef DO_HEAP2
    double comm1_max_key_old2;
#endif

    comm1_max_key_idx = communities[comm1].max_key_idx;
    assert(comm1_max_key_idx >= 0);
    comm1_max_key = communities[comm1].max_key;
    comm1_max_key_old = comm1_max_key;
#ifdef DO_HEAP2
    comm1_max_key_old2 = communities[comm1].adjcomm[comm1_max_key_idx].value[0];
#endif
    comm1_degree = communities[comm1].degree;
    comm1_adjcomm = communities[comm1].adjcomm;
    comm1_max_key_adjcomm_id = comm1_adjcomm[comm1_max_key_idx].comm_id;

#ifdef DEBUG_MB
    fprintf(stderr,"DEBUG_MB UpdateP2 %2d %2d\n",comm1, new_comm_id);
#endif

    heap  = maxheap->heap;
    index = maxheap->index;

    if (comm2 == comm1_max_key_adjcomm_id) {
	comm1_adjcomm[comm1_max_key_idx].value[KEY_CNM] = new_dq_val;
	comm1_adjcomm[comm1_max_key_idx].comm_id  = new_comm_id;
	if (new_lij.weighted)
	    comm1_adjcomm[comm1_max_key_idx].Xij      = new_lij.M;
	else
	    comm1_adjcomm[comm1_max_key_idx].lij      = new_lij.A;

#ifdef DEBUG_LIJ
	fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (update p2 comm2=comm1_max)\n",comm1,comm1_adjcomm[comm1_max_key_idx].comm_id, new_lij);
#endif

#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB bef UV %2d %2d p2 c1_max_key_idx: %2d \n",
		comm1, new_comm_id, comm1_max_key_idx);
#endif
	aggc_update_values(communities, comm1, comm1_max_key_idx, new_comm_id, L, maxheap->n - 1);

	memcpy(&new_adjcomm, comm1_adjcomm + comm1_max_key_idx, sizeof(aggc_adjcomm_t));

	if (comm1_adjcomm[comm1_max_key_idx].value[keytype] > comm1_max_key) {
	    comm1_max_key = comm1_adjcomm[comm1_max_key_idx].value[keytype];
	    communities[comm1].max_key = comm1_max_key;


	    /* move this adjacency to correct position */
	    if (new_comm_id > comm2) {
		for (i=comm1_max_key_idx; i<comm1_degree-1; i++) {
		    if (new_comm_id > comm1_adjcomm[i+1].comm_id) {
			/* swap */
			memcpy(comm1_adjcomm+i  , comm1_adjcomm+i+1, sizeof(aggc_adjcomm_t));
			memcpy(comm1_adjcomm+i+1, &new_adjcomm,      sizeof(aggc_adjcomm_t));
		    }
		    else {
			break;
		    }
		}
		communities[comm1].max_key_idx = i;
	    }
	    else {
		for (i=comm1_max_key_idx; i>0; i--) {
		    if (new_comm_id < comm1_adjcomm[i-1].comm_id) {
			/* swap */
			memcpy(comm1_adjcomm+i,   comm1_adjcomm+i-1, sizeof(aggc_adjcomm_t));
			memcpy(comm1_adjcomm+i-1, &new_adjcomm,      sizeof(aggc_adjcomm_t));
		    }
		    else {
			break;
		    }
		}
		communities[comm1].max_key_idx = i;
	    }

	}
	else {

	    /* move this adjacency to correct position */
	    if (new_comm_id > comm2) {
		for (i=comm1_max_key_idx; i<comm1_degree-1; i++) {
		    if (new_comm_id > comm1_adjcomm[i+1].comm_id) {
			/* swap */
			memcpy(comm1_adjcomm+i,   comm1_adjcomm+i+1, sizeof(aggc_adjcomm_t));
			memcpy(comm1_adjcomm+i+1, &new_adjcomm,      sizeof(aggc_adjcomm_t));
		    }
		    else {
			break;
		    }
		}
	    }
	    else {
		for (i=comm1_max_key_idx; i>0; i--) {
		    if (new_comm_id < comm1_adjcomm[i-1].comm_id) {
			/* swap */
			memcpy(comm1_adjcomm+i,   comm1_adjcomm+i-1, sizeof(aggc_adjcomm_t));
			memcpy(comm1_adjcomm+i-1, &new_adjcomm,      sizeof(aggc_adjcomm_t));
		    }
		    else {
			break;
		    }
		}
	    }
	    comm1_max_key     = NO_MAX_KEY;
	    comm1_max_key_idx = -1;
	    /* Rescan adjacencies to determine largest dq value */
	    for (i=0; i<comm1_degree; i++) {
		key_val = comm1_adjcomm[i].value[keytype];
		if (key_val > comm1_max_key) {
		    comm1_max_key_idx = i;
		    comm1_max_key     = key_val;
		}
	    }
	    communities[comm1].max_key     = comm1_max_key;
	    communities[comm1].max_key_idx = comm1_max_key_idx;
	}

	/* Sift up/down value in max heap */
	heap_pos = index[comm1];

	/* heap_pos < 0 is only possible during seed-set growth, not
	   the global alg. */
	if (heap_pos >= 0) {
	    heap[heap_pos].val = comm1_max_key;

#ifdef DO_HEAP2
	    if (maxheap->use2)
		heap[heap_pos].val2 = comm1_adjcomm[comm1_max_key_idx].value[0];
#endif

	    if ( (fpGT(comm1_max_key,comm1_max_key_old))
#ifdef DO_HEAP2
		 || ((maxheap->use2)&&(fpEQ(comm1_max_key,comm1_max_key_old))&&(fpGT(heap[heap_pos].val2,comm1_max_key_old2)))
#endif
		 )
		aggc_maxheap_sift_up(maxheap, heap_pos);
	    else
		aggc_maxheap_sift_down(maxheap, heap_pos);

#ifdef DEBUG_MAXHEAP
	    fprintf(stderr,"DEBUG_MAXHEAP Checking heap: p2 A\n");
	    aggc_maxheap_check(maxheap);
#endif
	}
    }
    else {

	/* Locate this adjacency */
	pos = aggc_adjcomm_pos(comm1_adjcomm, comm2, comm1_degree);
	comm1_adjcomm[pos].value[KEY_CNM] = new_dq_val;
	comm1_adjcomm[pos].comm_id  = new_comm_id;
	if (new_lij.weighted)
	    comm1_adjcomm[pos].Xij      = new_lij.M;
	else
	    comm1_adjcomm[pos].lij      = new_lij.A;

#ifdef DEBUG_LIJ
	fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (update p2-else)\n",comm1,comm1_adjcomm[pos].comm_id, new_lij);
#endif

	memcpy(&new_adjcomm, comm1_adjcomm + pos, sizeof(aggc_adjcomm_t));

	/* move this adjacency to correct position */
	if (new_comm_id > comm2) {
	    for (i=pos; i<comm1_degree-1; i++) {
		if (new_comm_id > comm1_adjcomm[i+1].comm_id) {
		    /* swap */
		    memcpy(comm1_adjcomm+i,   comm1_adjcomm+i+1, sizeof(aggc_adjcomm_t));
		    memcpy(comm1_adjcomm+i+1, &new_adjcomm,      sizeof(aggc_adjcomm_t));
		}
		else {
		    break;
		}
	    }
	}
	else {
	    for (i=pos; i>0; i--) {
		if (new_comm_id < comm1_adjcomm[i-1].comm_id) {
		    /* swap */
		    memcpy(comm1_adjcomm+i,   comm1_adjcomm+i-1, sizeof(aggc_adjcomm_t));
		    memcpy(comm1_adjcomm+i-1, &new_adjcomm,      sizeof(aggc_adjcomm_t));
		}
		else {
		    break;
		}
	    }
	}

#if 1


#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB bef UV %2d %2d p2-else pos: %2d i: %2d\n",
		comm1, new_comm_id, pos, i);
#endif
	aggc_update_values(communities, comm1, i, new_comm_id, L, maxheap->n - 1);
#endif

	/* Update its key value */
	if (comm1_adjcomm[pos].value[keytype] > comm1_max_key) {
	    /* Update global adj value */
	    communities[comm1].max_key     = comm1_adjcomm[pos].value[keytype];
	    communities[comm1].max_key_idx = i;

	    /* Sift up/down value in max heap */
	    heap_pos = index[comm1];

	    /* heap_pos < 0 is only possible during seed-set growth,
	       not the global alg. */
	    if (heap_pos >= 0) {
		heap[heap_pos].val = communities[comm1].max_key;
#ifdef DO_HEAP2
		if (maxheap->use2)
		    heap[heap_pos].val2 = communities[comm1].adjcomm[communities[comm1].max_key_idx].value[0];
#endif
		aggc_maxheap_sift_up(maxheap, heap_pos);
#ifdef DEBUG_MAXHEAP
		fprintf(stderr,"DEBUG_MAXHEAP Checking heap: p2 B\n");
		aggc_maxheap_check(maxheap);
#endif
	    }
	}
	else {
	    /* update local comm1_max_key_idx */
	    comm1_max_key_idx = communities[comm1].max_key_idx;
	    if (new_comm_id > comm2) {
		if ((comm1_max_key_idx > pos) && (comm1_max_key_idx <= i)) {
		    communities[comm1].max_key_idx = comm1_max_key_idx - 1;
		}
	    }
	    else {
		if ((comm1_max_key_idx >= i) && (comm1_max_key_idx < pos)) {
		    communities[comm1].max_key_idx = comm1_max_key_idx + 1;
		}
	    }
	}
    }

    return;
}


static void
aggc_remove_adjcomm(aggc_comm_t *communities, aggc_maxheap_t *maxheap,
		    attr_id_t comm1, attr_id_t comm2)
{

    attr_id_t i, comm_id, max_key_idx, max_key_adjcomm_id, pos, degree;
    aggc_adjcomm_t *adjcomm;
    aggc_heap_t *heap;
    double max_key, max_key_old, key_val;
    attr_id_t *index;
    attr_id_t heap_pos;
#ifdef DO_HEAP2
    double max_key_old2;
#endif

    index   = maxheap->index;
    heap    = maxheap->heap;

    max_key_old = communities[comm1].max_key;
#ifdef DO_HEAP2
    max_key_old2 = communities[comm1].adjcomm[communities[comm1].max_key_idx].value[0];
#endif
    max_key_idx = communities[comm1].max_key_idx;

    if (max_key_idx < 0) {
#if 0
	fprintf(stderr, "comm %d, max_key_idx %d, degree %d, dqi %d\n",
		comm1, max_key_idx, index[comm1]);
#endif
	assert(max_key_idx >= 0);
    }
    max_key_adjcomm_id = communities[comm1].adjcomm[max_key_idx].comm_id;
    degree = communities[comm1].degree;
    adjcomm = communities[comm1].adjcomm;

    if (comm2 == max_key_adjcomm_id) {

	/* we need to rescan to find new max_key for comm1 */
	max_key     = NO_MAX_KEY;
	max_key_idx = -1;
	communities[comm1].max_key     = NO_MAX_KEY;
	communities[comm1].max_key_idx = -1;

	i = 0;
	while (i < degree) {
	    comm_id = adjcomm[i].comm_id;
	    key_val = adjcomm[i].value[keytype];
	    if (comm_id == comm2)
		break;
	    if (key_val > max_key) {
		max_key = key_val;
		max_key_idx = i;
	    }
	    i++;
	}

	pos = i;
	for (i=pos; i<degree-1; i++) {
	    memcpy(adjcomm+i, adjcomm+i+1, sizeof(aggc_adjcomm_t));
	    if (adjcomm[i].value[keytype] > max_key) {
		max_key = adjcomm[i].value[keytype];
		max_key_idx = i;
	    }
	}

	communities[comm1].max_key = max_key;
	communities[comm1].max_key_idx = max_key_idx;

	/* Sift down value in max-heap */
	heap_pos = index[comm1];

	/* heap_pos < 0 is only possible during seed-set growth, not
	   the global alg. */
	if (heap_pos >= 0) {
	    heap[heap_pos].val = max_key;
#ifdef DO_HEAP2
	    if (maxheap->use2)
		heap[heap_pos].val2 = communities[comm1].adjcomm[communities[comm1].max_key_idx].value[0];
#endif

	    if ( (fpGT(max_key,max_key_old))
#ifdef DO_HEAP2
		 || ((maxheap->use2)&&(fpEQ(max_key,max_key_old))&&(fpGT(heap[heap_pos].val2,max_key_old2)))
#endif
		 )
		aggc_maxheap_sift_up(maxheap, heap_pos);
	    else
		aggc_maxheap_sift_down(maxheap, heap_pos);

#ifdef DEBUG_MAXHEAP
	    fprintf(stderr,"DEBUG_MAXHEAP Checking heap: remove adjcomm A\n");
	    aggc_maxheap_check(maxheap);
#endif
	}
    }
    else {

	/* do a binary search to locate the pair in adjcomm, delete it,
	 * and update the communities ordering */
	pos = aggc_adjcomm_pos(adjcomm, comm2, degree);
	if (max_key_idx > pos)
	    communities[comm1].max_key_idx--;

	for (i=pos; i<degree-1; i++) {
	    memcpy(adjcomm+i, adjcomm+i+1, sizeof(aggc_adjcomm_t));
	}
    }

    communities[comm1].degree--;
    return;
}

void aggc_compute_dqnorm_stats(aggc_comm_t *communities, aggc_maxheap_t *maxheap,
	double *dqnorm_mean, double *dqnorm_stddev, attr_id_t *edgeCnt) {
    attr_id_t i, j, cnt;
    double x, sumx, sumx2, mean;

    cnt = 0;
    sumx  = 0.0;
    sumx2 = 0.0;
    for (i=0; i<maxheap->n; i++) {
	for (j=0 ; j<communities[i].degree ; j++) {
	    x = communities[i].adjcomm[j].value[keytype];
#ifdef DEBUG_STATS
	    fprintf(stderr,"DEBUG_STATS x(%5d, %5d): [%9.6f  %9.6f]\n",i,j,communities[i].adjcomm[j].value[KEY_CNM],x);
#endif
	    sumx += x;
	    sumx2 += x*x;
	    cnt++;
	}
    }
    *dqnorm_mean   = mean = sumx/(double)cnt;
    *dqnorm_stddev = sqrt(sumx2/(double)cnt - (mean * mean));
#ifdef DEBUG_STATS
    fprintf(stderr,"DEBUG_STATS sumx2/num_comm: %9.6f  mean*mean: %9.6f\n",sumx2/(double)cnt, mean*mean);
    fprintf(stderr,"DEBUG_STATS mean: %9.6f  stddev: %9.6f   cnt: %5d\n",mean, *dqnorm_stddev, cnt);
#endif
    *edgeCnt = cnt;

    return;
}

double
aggc_merge_communities(aggc_comm_t *communities,
		       aggc_maxheap_t *maxheap,
		       aggc_adjcomm_t *adj_buffer,
		       attr_id_t num_communities,
		       weight_t *L)
{
    attr_id_t HIGH_VID;
    attr_id_t p1, p2, comm1, comm2, from, to, i, j, from_deg, to_deg;
    aggc_adjcomm_t *from_adj, *to_adj;
    aggc_heap_t *heap;
    double max_dq, max_key, from_a, to_a, max_key_old, dqnorm_mean, dqnorm_stddev, to_max_key;
    double max_edge_2Xij;
    attr_id_t ins_pos, max_key_idx, to_max_key_idx, curr_i_idx, heap_pos, *index, edgeCnt;
    double tmp, MKH_mean, MKH_stddev;
    weight_t updated_weight;
#ifdef DO_HEAP2
    double max_key_old2;
#endif

    index   = maxheap->index;
    heap    = maxheap->heap;

    updated_weight.weighted = L->weighted;

    HIGH_VID = 1<<30;
    comm1 = heap[0].comm_id;
    max_key = heap[0].val;
    max_key_idx = communities[comm1].max_key_idx;

    if (max_key != communities[comm1].max_key) {
#if 1
	fprintf(stderr,"ERROR: max_key != communities[comm1].max_key, ");
	fprintf(stderr,"ERROR: merge: comm1: %5d  max_key: %f  comm[comm1].max_key: %f  max_key_idx: %d\n",
		comm1, max_key, communities[comm1].max_key, max_key_idx);
#endif
	return 0;
    }


#ifdef DEBUG_MAXHEAP
    fprintf(stderr,"DEBUG_MAXHEAP Checking heap: start of merge comm\n");
    aggc_maxheap_check(maxheap);
#endif

    assert(max_key == communities[comm1].max_key);

#ifdef DEBUG
    fprintf(stderr,"DEBUG merge: max_key: %f  comm[comm1].max_key: %f  idx: %d\n",
	    max_key, communities[comm1].max_key, max_key_idx);
    /*  assert(max_key == communities[comm1].max_key); */
    fprintf(stderr,"DEBUG merge: communities[comm1].adjcomm[max_key_idx].value[KEY_CNM]: %f\n",communities[comm1].adjcomm[max_key_idx].value[KEY_CNM]);
#endif

    max_dq = communities[comm1].adjcomm[max_key_idx].value[KEY_CNM];

    /* Stopping criteron */
    switch (keytype) {

	case KEY_CNM:
#ifdef DEBUG_STATS
	    fprintf(stderr,"DEBUG_STATS stopping check  max_key: %9.6f\n", max_key);
#endif
	    if (max_dq < 0.0) {
#ifdef VERBOSE
		fprintf(stderr,"Stopping: max_dq: %f < 0.\n", max_dq);
#endif
		return max_dq;
	    }
	    break;

	case KEY_MB:

#if 0
	    aggc_compute_dqnorm_stats(communities, maxheap, &dqnorm_mean, &dqnorm_stddev, &edgeCnt);
	    /* if the max_dqnorm is not statistically significant */
#ifdef DEBUG_STATS
	    fprintf(stderr,"DEBUG_STATS dqnorm (%9.6f, %9.6f)  max_key: %9.6f\n",
		    dqnorm_mean, dqnorm_stddev, max_key);
	    fprintf(stderr,"DEBUG_STATS stop e/L: %9.6f  mkey: %9.6f stat: %9.6f  mdq: %9.6f\n",
		    (double)edgeCnt/(double)(2*L->A),
		    max_key,
		    dqnorm_mean - NSTDDEV*dqnorm_stddev,
		    max_dq);
#endif
	    if ( (max_key < dqnorm_mean - NSTDDEV*dqnorm_stddev) || (max_dq < 0) ) {
		return max_dq;
	    }
#else
	    if (max_dq < 0) {
#ifdef VERBOSE
		fprintf(stderr,"Stopping: max_dq: %f < 0.\n", max_dq);
#endif
		return max_dq;
	    }

	    if (max_key < INFTY) {
		tmp = MaxKeyHistory[MaxKeyHistoryIdx];
		MaxKeyHistory[MaxKeyHistoryIdx] = max_key;
		MaxKeyHistoryIdx = (MaxKeyHistoryIdx + 1) % MAXHISTORY;

		if (MaxKeyHistoryFilled) {
		    sumMKH  -= tmp;
		    sumMKH2 -= tmp*tmp;
		}
		sumMKH  += max_key;
		sumMKH2 += max_key*max_key;

		if ((!MaxKeyHistoryFilled)&&(MaxKeyHistoryIdx==0))
		    MaxKeyHistoryFilled = 1;

		MKH_mean = sumMKH/(double)MAXHISTORY;
		MKH_stddev = sqrt(sumMKH2/(double)MAXHISTORY - (MKH_mean * MKH_mean));
#if 0
		fprintf(stderr,"stop e/L: %9.6f  mkey: %9.6f stat: %9.6f  mdq: %9.6f\n",
			(double)edgeCnt/(double)(2*L->A),
			max_key,
			MKH_mean - NSTDDEV*MKH_stddev,
			max_dq);
#endif

#ifdef DO_STOPPING_LIJ
		if (MaxKeyHistoryFilled && (max_key < MKH_mean - NSTDDEV*MKH_stddev)) {
#if 0 || defined(VERBOSE)
		    fprintf(stderr,"Stopping: max_key: %f < mean - %g s.d.'s\n", max_key, NSTDDEV);
#endif
		    return max_dq;
		}
#endif
	    }

#endif
	    break;

	case KEY_RAT:
	case KEY_MBRAT:
	case KEY_LIN:

#if 0
	    if (max_dq < 0.0) {
#ifdef VERBOSE
		fprintf(stderr,"Stopping: max_dq: %f < 0.\n", max_dq);
#endif
		return max_dq;
	    }
#endif

#if 0
	    if (index[comm1]<0) {
#ifdef VERBOSE
		fprintf(stderr,"Stopping: index[comm1] < 0.\n");
#endif
		return max_dq;
	    }
#endif


#if 0
	    if (index[communities[comm1].adjcomm[max_key_idx].comm_id]<0) {
#ifdef VERBOSE
		fprintf(stderr,"Stopping: index[comm2] < 0.\n");
#endif
		return max_dq;
	    }
#endif

#if 0
	    {
		attr_id_t c2;
		double edgeWt, c1density, c2density;

		c2 = communities[comm1].adjcomm[max_key_idx].comm_id;
		edgeWt = communities[comm1].adjcomm[max_key_idx].Xij;

		c1density = communities[comm1].comm_wt / (double)communities[comm1].comm_size;
		c2density = communities[c2   ].comm_wt / (double)communities[c2   ].comm_size;

#ifdef DEBUG
		fprintf(stderr,"CHECK c1: %5d c2: %5d wt: %12f c1w: %12f c1a: %12f c1d: %12f c2w: %12f c2a: %12f c2d: %12f ratio: %12f newdens: %12f\n",
			comm1, c2, edgeWt,
			communities[comm1].comm_wt,
			communities[comm1].a_lin * L->origM,
			c1density,
			communities[c2].comm_wt,
			communities[c2].a_lin * L->origM,
			c2density,
			(c1density > c2density ? c1density/c2density : c2density/c1density),
			(communities[comm1].comm_wt + 2.0*edgeWt + communities[c2].comm_wt)/(double)(communities[comm1].comm_size+communities[c2].comm_size));
#endif

		if ( (communities[comm1].comm_wt > communities[comm1].a_lin * L->origM) && (communities[c2].comm_wt > communities[c2].a_lin * L->origM)) {
#ifdef VERBOSE
		    fprintf(stderr,"Stopping: communities both have more internal than external weight\n");
#endif
		    return max_dq;
		}

		if ((communities[comm1].comm_size > 4) && (communities[c2].comm_size > 4) && (c1density > 0.0 && c2density > 0.0)) {

#ifdef DEBUG
		    fprintf(stderr,"CHECKSTOP? wt: %12f val: %12f c1w: %12f c1d: %12f c2w: %12f c2d: %12f ratio: %12f newdens: %12f <--\n",
			    edgeWt, max_key,
			    communities[comm1].comm_wt,
			    c1density,
			    communities[c2].comm_wt,
			    c2density,
			    (c1density > c2density ? c1density/c2density : c2density/c1density),
			    (communities[comm1].comm_wt + 2.0*edgeWt + communities[c2].comm_wt)/(double)(communities[comm1].comm_size+communities[c2].comm_size));
#endif

#if 0
		    if ((edgeWt < communities[comm1].comm_wt) && (edgeWt < communities[c2].comm_wt)) {
#ifdef VERBOSE
			fprintf(stderr,"Stopping: edgeWt too small: %f\n", edgeWt);
#endif
			return max_dq;
		    }
#endif

#if 0
		    if (edgeWt <  (c1density > c2density ? c1density/c2density : c2density/c1density) )
			fprintf(stderr,"STOP CHECK: weight less than than ratio\n");


		    fprintf(stderr,"STOP CHECK: Xij (%f) not significant to comm1 size: %d wt: %f dens: %f and comm2 size: %d wt: %f dens: %f.\n",
			    edgeWt,
			    communities[comm1].comm_size,
			    communities[comm1].comm_wt,
			    c1density,
			    communities[c2].comm_size,
			    communities[c2].comm_wt,
			    c2density);

#endif
		}


#if 0
		if ((edgeWt <= c1density / 4.0) &&  (edgeWt <= c2density / 4.0)) {
#ifdef VERBOSE
		    fprintf(stderr,"Stopping: Xij (%f) not significant to comm1 size: %d and comm2 size: %d.\n",
			    edgeWt,
			    communities[comm1].comm_size,
			    communities[c2   ].comm_size);
#endif
		    return max_dq;
		}
#endif

	    }
#endif

#if 1
	    if (max_key < 0.0) {
#ifdef VERBOSE
		fprintf(stderr,"Stopping: max_key: %f < 0.\n", max_key);
#endif
		return max_dq;
	    }
#endif

	    break;

	default:
	    fprintf(stderr,"No stopping criteron\n");
	    exit(-1);
    }

    comm2 = communities[comm1].adjcomm[max_key_idx].comm_id;

    if (comm1 != communities[comm2].adjcomm[communities[comm2].max_key_idx].comm_id) {
	/* comm2 has an adjacency comm3 which has the same dq value
	   as the pair comm1 -> comm2 */
	/* Search to locate comm1 and modify comm2's max_key_idx */
	ins_pos = aggc_adjcomm_pos(communities[comm2].adjcomm, comm1, communities[comm2].degree);
	communities[comm2].max_key_idx = ins_pos;
    }

    /* If this is a directed merge for seed-set expansion, index[comm2] < 0. */
    if (index[comm2] >= 0 && communities[comm1].degree < communities[comm2].degree) {
	from = comm1;
	to = comm2;
    }
    else {
	from = comm2;
	to = comm1;
    }

    from_deg = communities[from].degree;
    from_adj = communities[from].adjcomm;
    from_a   = communities[from].a;
    to_deg   = communities[to].degree;
    to_adj   = communities[to].adjcomm;
    to_a     = communities[to].a;

#ifdef DEBUG_MAXKEY
    /* Make sure we haven't violated to's maxkey! */
    fprintf(stderr,"DEBUG Maxkey check for _to_ BEFORE.\n");
    aggc_maxkey_comm_check(communities, to);
    fprintf(stderr,"DEBUG maxkey to BEFORE to: %d  degree: %d  max_key: %f  max_key_idx: %d\n",to,communities[to].degree,
	    communities[to].max_key, communities[to].max_key_idx);
    for (i=0 ; i<communities[to].degree; i++)
	fprintf(stderr,"DEBUG maxkey to BEFORE adjcomm[%2d].value[keytype]: %f  (%2d comm_id)\n",
		i, communities[to].adjcomm[i].value[keytype], communities[to].adjcomm[i].comm_id);
#endif


    to_max_key     = NO_MAX_KEY;
    to_max_key_idx = -1;
    curr_i_idx     = 0;

#ifdef DEBUG
    fprintf(stderr, "DEBUG Merging communities %12d -> %12d   max_dq %9.6f\n", from, to, max_dq);
#endif
#ifdef VERBOSE
    fprintf(stdout, "(n: %8d) Merging comm's %8d (d: %8d, s: %8d) -> %8d (d: %8d, s: %8d), mod-inc: %12.9f\n",
	    maxheap->n,
	    from, communities[from].degree, communities[from].comm_size,
	    to,communities[to].degree, communities[to].comm_size,
	    max_dq);
#endif

#if 0
    if ((keytype==KEY_LIN)&&(max_dq < 0.0)) {
	fprintf(stdout,"POTENTIAL T NODE: %12d\n",from);
    }
#endif

    if (keytype==KEY_LIN) {
	assert (L->weighted);
	max_edge_2Xij = 2.0 * communities[comm1].adjcomm[max_key_idx].Xij;
	L->M -= max_edge_2Xij;
#ifdef DEBUG_LIN
	fprintf(stderr,"DEBUG_LIN reducing L->M by 2*%f, new L->M: %f\n",
		communities[comm1].adjcomm[max_key_idx].Xij, L->M);
#endif
    }

    i = 0;
    j = 0;
    while (i<to_deg && j<from_deg) {

	p1 = to_adj[i].comm_id;
	p2 = from_adj[j].comm_id;

#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB Merging adjacencies to: %2d i: %2d from: %2d j: %2d p1: %2d p2: %2d\n",
		to, i, from, j, p1, p2);
#endif

	if (p1 == from) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (p1==from) to: %2d i: %2d from: %2d j: %2d p1: %2d p2: %2d\n",
		    to, i, from, j, p1, p2);
#endif
	    /* This is "to" pointing to "from" */
	    /* This community needs to be deleted */
	    to_adj[i].comm_id = HIGH_VID;
	    i++;
	    continue;
	}

	if (p2 == to) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (p2==to) to: %2d i: %2d from: %2d j: %2d p1: %2d p2: %2d\n",
		    to, i, from, j, p1, p2);
#endif
	    /* This is "from" pointing to "to" */
	    j++;
	    continue;
	}

	if (p1 < p2) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (p1<p2) to: %2d i: %2d from: %2d j: %2d p1: %2d p2: %2d\n",
		    to, i, from, j, p1, p2);
#endif
	    /* We have a community pair that is in "to" but not in "from" */
	    /* Update its dq value */
	    to_adj[i].value[KEY_CNM] -= 2*from_a*communities[p1].a;

#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB bef UV %2d %2d mc p1<p2 to: %2d i: %2d  p1: %2d\n",
		    to, p1, to, i, p1);
#endif
	    aggc_update_values(communities, to, i, p1, *L, maxheap->n - 1);

	    if (to_adj[i].value[keytype] > to_max_key) {
		to_max_key = to_adj[i].value[keytype];
		to_max_key_idx = curr_i_idx;
	    }

	    /* Update dq value of adjcomm (p1 -> to) */
	    if (updated_weight.weighted)
		updated_weight.M = to_adj[i].Xij;
	    else
		updated_weight.A = to_adj[i].lij;

	    aggc_update_dq_p1(communities, maxheap, p1, to, *L, to_adj[i].value[KEY_CNM], updated_weight);
	    i++;
	    curr_i_idx++;
	    continue;
	}

	if (p1 == p2) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (p1==p2) to: %2d i: %2d from: %2d j: %2d p1: %2d p2: %2d\n",
		    to, i, from, j, p1, p2);
#endif

	    /* We have a community that is adjacent to both "from" and "to" */
	    to_adj[i].value[KEY_CNM] += from_adj[j].value[KEY_CNM];
	    if (updated_weight.weighted)
		to_adj[i].Xij += from_adj[j].Xij;
	    else
		to_adj[i].lij += from_adj[j].lij;
#ifdef DEBUG_LIJ
	    fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (merge p1==p2)\n",to,to_adj[i].comm_id, to_adj[i].lij);
#endif


#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB bef UV %2d %2d mc p1==p2 to: %2d i: %2d  p1: %2d\n",
		    to, p1, to, i, p1);
#endif
	    aggc_update_values(communities, to, i, p1, *L, maxheap->n - 1);

	    if (to_adj[i].value[keytype] > to_max_key) {
		to_max_key     = to_adj[i].value[keytype];
		to_max_key_idx = curr_i_idx;
	    }

	    if (updated_weight.weighted)
		updated_weight.M = to_adj[i].Xij;
	    else
		updated_weight.A = to_adj[i].lij;

	    aggc_update_dq_p1(communities, maxheap, p1, to, *L, to_adj[i].value[KEY_CNM], updated_weight);

	    aggc_remove_adjcomm(communities, maxheap, p2, from);

#ifdef DEBUG_MAXKEY
	    /* Make sure we haven't violated p1's maxkey! */
	    fprintf(stderr,"DEBUG_MAXKEY Maxkey check for p1.\n");
	    aggc_maxkey_comm_check(communities, p1);
#endif

	    i++;
	    j++;
	    curr_i_idx++;
	    continue;
	}

	if (p1 > p2) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (p1>p2) to: %2d i: %2d from: %2d j: %2d p1: %2d p2: %2d\n",
		    to, i, from, j, p1, p2);
#endif
	    /* We have a community pair that is in "from" but not in "to" */
	    if (communities[to].degree == communities[to].max_size) {
		aggc_adjcomm_resize(communities, to, from_deg);
		to_adj = communities[to].adjcomm;
	    }
	    ins_pos = communities[to].degree;
	    communities[to].degree++;
	    to_adj[ins_pos].comm_id = p2;
	    to_adj[ins_pos].value[KEY_CNM] = from_adj[j].value[KEY_CNM] - 2*to_a*communities[p2].a;
	    if (updated_weight.weighted)
		to_adj[ins_pos].Xij = from_adj[j].Xij;
	    else
		to_adj[ins_pos].lij = from_adj[j].lij;
#ifdef DEBUG_LIJ
	    fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (merge p1>p2)\n",to,to_adj[ins_pos].comm_id, to_adj[ins_pos].lij);
#endif

	    /* compute other values of to_adj[ins_pos] */
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB bef UV %2d %2d mc p1>p2 to: %2d ins_pos: %2d  p2: %2d\n",
		    to, p2, to, ins_pos, p2);
#endif
	    aggc_update_values(communities, to, ins_pos, p2, *L, maxheap->n - 1);

	    if (to_adj[ins_pos].value[keytype] > to_max_key) {
		to_max_key = to_adj[ins_pos].value[keytype];
		to_max_key_idx = curr_i_idx;
	    }

	    /* Update (p2 -> from) */
	    if (updated_weight.weighted)
		updated_weight.M = to_adj[i].Xij;
	    else
		updated_weight.A = to_adj[ins_pos].lij;

	    aggc_update_dq_p2(communities, maxheap, p2, from, to, *L, to_adj[ins_pos].value[KEY_CNM], updated_weight);


#ifdef DEBUG_MAXKEY
	    /* Make sure we haven't violated p2's maxkey! */
	    fprintf(stderr,"DEBUG_MAXKEY Maxkey check for p2.\n");
	    aggc_maxkey_comm_check(communities, p2);
#endif

	    j++;
	    curr_i_idx++;
	    continue;
	}

    }

    while (i < to_deg) {
	p1 = to_adj[i].comm_id;
#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB Merging (i<to_deg) to: %2d i: %2d from: %2d j: %2d p1: %2d \n",
		to, i, from, j, p1);
#endif
	if (p1 == from) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (i<to_deg),(p1==from) to: %2d i: %2d from: %2d j: %2d p1: %2d \n",
		    to, i, from, j, p1);
#endif
	    /* This community needs to be deleted */
	    to_adj[i].comm_id = HIGH_VID;
	    i++;
	    continue;
	}
	to_adj[i].value[KEY_CNM] -= 2*from_a*communities[p1].a;

	/* update other values of to_adj[i] */
#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB bef UV %2d %2d while-i to: %2d i: %2d  p1: %2d\n",
		to, p1, to, i, p1);
#endif
	aggc_update_values(communities, to, i, p1, *L, maxheap->n - 1);

#ifdef DEBUG_MAXKEY
	/* Make sure we haven't violated p1's maxkey! */
	fprintf(stderr,"DEBUG_MAXKEY Maxkey check for p1.\n");
	aggc_maxkey_comm_check(communities, p1);
#endif

	if (to_adj[i].value[keytype] > to_max_key) {
	    to_max_key = to_adj[i].value[keytype];
	    to_max_key_idx = curr_i_idx;
	}

	/* Update dq value of adjcomm (p1 -> to) */
	if (updated_weight.weighted)
	    updated_weight.M = to_adj[i].Xij;
	else
	    updated_weight.A = to_adj[i].lij;

	aggc_update_dq_p1(communities, maxheap, p1, to, *L, to_adj[i].value[KEY_CNM], updated_weight);
	i++;
	curr_i_idx++;
    }

    while (j < from_deg) {
	p2 = from_adj[j].comm_id;
#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB Merging (j<from_deg) to: %2d i: %2d from: %2d j: %2d p2: %2d\n",
		to, i, from, j, p2);
#endif
	if (p2 == to) {
#ifdef DEBUG_MB
	    fprintf(stderr,"DEBUG_MB Merging (j<from_deg)(p2==to) to: %2d i: %2d from: %2d j: %2d p2: %2d\n",
		    to, i, from, j, p2);
#endif
	    j++;
	    continue;
	}
	if (communities[to].degree == communities[to].max_size) {
	    aggc_adjcomm_resize(communities, to, from_deg);
	    to_adj = communities[to].adjcomm;
	}
	ins_pos = communities[to].degree;
	communities[to].degree++;
	to_adj[ins_pos].comm_id = p2;
	to_adj[ins_pos].value[KEY_CNM] = from_adj[j].value[KEY_CNM] - 2*to_a*communities[p2].a;
	if (updated_weight.weighted)
	    to_adj[ins_pos].Xij = from_adj[j].Xij;
	else
	    to_adj[ins_pos].lij = from_adj[j].lij;
#ifdef DEBUG_LIJ
	fprintf(stderr,"DEBUG_LIJ update %2d %2d <-- %2d (merge j<from_deg)\n",to,to_adj[ins_pos].comm_id, to_adj[ins_pos].lij);
#endif

	/* update other values of to_adj[ins_pos] */
#ifdef DEBUG_MB
	fprintf(stderr,"DEBUG_MB bef UV %2d %2d while-j to: %2d ins_pos: %2d  p2: %2d\n",
		to, p2, to, ins_pos, p2);
#endif
	aggc_update_values(communities, to, ins_pos, p2, *L, maxheap->n - 1);

#ifdef DEBUG_MAXKEY
	/* Make sure we haven't violated p2's maxkey! */
	fprintf(stderr,"DEBUG_MAXKEY Maxkey check for p2.\n");
	aggc_maxkey_comm_check(communities, p2);
#endif

	if (to_adj[ins_pos].value[keytype] > to_max_key) {
	    to_max_key = to_adj[ins_pos].value[keytype];
	    to_max_key_idx = curr_i_idx;
	}

	/* Update (p2 -> from) */
	if (updated_weight.weighted)
	    updated_weight.M = to_adj[ins_pos].Xij;
	else
	    updated_weight.A = to_adj[ins_pos].lij;

	aggc_update_dq_p2(communities, maxheap, p2, from, to, *L, to_adj[ins_pos].value[KEY_CNM], updated_weight);
	j++;
	curr_i_idx++;
    }

    assert(i == to_deg);
    assert(j == from_deg);

#ifdef DEBUGXX
    for (i=0; i<communities[to].degree; i++) {
	fprintf(stderr, "DEBUGXX Merging TO (%3d): [%d %ld]\n", to, communities[to].adjcomm[i].comm_id,
		communities[to].adjcomm[i].value[KEY_CNM]);
    }
#endif

    qsort(communities[to].adjcomm, communities[to].degree, sizeof(aggc_adjcomm_t), comp_comm);

    communities[to].a += communities[from].a;

    if (keytype==KEY_LIN) {
#ifdef DO_ALIN
	communities[to].a_lin += communities[from].a_lin;
	communities[to].a_lin -= max_edge_2Xij/L->origM;
#ifdef DO_COMM_WT
	communities[to].comm_wt += communities[from].comm_wt;
	communities[to].comm_wt += max_edge_2Xij;
#endif
#else
	communities[to].a -= max_edge_2Xij/L->origM;
#endif
    }

    communities[to].degree -= 1;
#if 0
    communities[to].comm_size += 1;
#else
    communities[to].comm_size += communities[from].comm_size;
#endif

    communities[from].parent_id = to;

    if (communities[from].resized == 1) {
	free(communities[from].adjcomm);
	communities[from].resized = 0;
    }

#ifdef DEBUGXX
    heap_pos = index[from];
    if (max_key_old != heap[heap_pos].value[KEY_CNM]) {
	fprintf(stderr, "DEBUGXX %lf %lf\n", max_key_old, heap[heap_pos].val);
    }
#endif

    aggc_maxheap_remove(maxheap, from);


#ifdef DEBUG_MAXKEY
    fprintf(stderr,"DEBUG_MAXKEY old maxkey: %f  old max_key_idx: %2d\n",
	    communities[to].max_key,
	    communities[to].max_key_idx);
#endif

#ifdef CULL_DEGREE_ONE
    if (keytype != KEY_CNM) {
	if ((communities[to].degree == 1) && (communities[to].comm_size <= CULL_THRESHOLD) ) {
#ifdef VERBOSE
	    fprintf(stderr,"Culling community %12d (size %12d):\n",to,communities[to].comm_size);
#endif
	    communities[to].adjcomm[0].value[keytype] = INFTY;
	    to_max_key = INFTY;
	    to_max_key_idx = 0;
	}
    }
#endif

    communities[to].max_key     = to_max_key;
    communities[to].max_key_idx = to_max_key_idx;


#ifdef DEBUG_MAXKEY
    /* Make sure we haven't violated to's maxkey! */
    fprintf(stderr,"DEBUG_MAXKEY Maxkey check for _to_ AFTER.\n");
    aggc_maxkey_comm_check(communities, to);
    fprintf(stderr,"DEBUG_MAXKEY maxkey to AFTER to: %d  degree: %d  max_key: %f  max_key_idx: %d\n",to,communities[to].degree,
	    communities[to].max_key, communities[to].max_key_idx);
    for (i=0 ; i<communities[to].degree; i++)
	fprintf(stderr,"DEBUG_MAXKEY maxkey to AFTER adjcomm[%2d].value[keytype]: %f  (%2d comm_id)\n",
		i, communities[to].adjcomm[i].value[keytype], communities[to].adjcomm[i].comm_id);
#endif

    heap_pos = index[to];
    assert(heap_pos >= 0);

    max_key_old =  heap[heap_pos].val;
#ifdef DO_HEAP2
    max_key_old2 = heap[heap_pos].val2;
#endif
    heap[heap_pos].val = to_max_key;
#ifdef DO_HEAP2
    heap[heap_pos].val2 = communities[to].adjcomm[to_max_key_idx].value[0];
#endif
    if ( (fpGT(to_max_key,max_key_old))
#ifdef DO_HEAP2
	    || ((maxheap->use2)&&(fpEQ(to_max_key,max_key_old))&&(fpGT(heap[heap_pos].val2,max_key_old2)))
#endif
       )
	aggc_maxheap_sift_up(maxheap, heap_pos);
    else
	aggc_maxheap_sift_down(maxheap, heap_pos);

#ifdef DEBUG_MAXHEAP
    fprintf(stderr,"DEBUG_MAXHEAP Checking heap: remove adjcomm B\n");
    aggc_maxheap_check(maxheap);
#endif

    return max_dq;
}



attr_id_t
aggc_compute_membership(aggc_comm_t *communities, attr_id_t *membership, attr_id_t n)
{
    attr_id_t final_num_communities, i, parent_id;

    /* get final community membership information */
    final_num_communities = 0;
    for (i=0; i<n; i++)
	if (communities[i].parent_id == -1)
	    membership[i] = final_num_communities++;

    for (i=0; i<n; i++) {
	if (communities[i].parent_id != -1) {
	    parent_id = communities[i].parent_id;
	    while (communities[parent_id].parent_id != -1)
		parent_id = communities[parent_id].parent_id;
	    membership[i] = membership[parent_id];
	}
    }

    return final_num_communities;
}




/* Greedy community detection algorithms that optimize modularity. These are
 * agglomerative approaches, i.e., we start by assuming a partitioning of
 * singleton communities, and then repeatedly merge communities with the
 * goal of maximizing modularity. There are several heuristics for this, and we
 * use a similar representation of communities for all of them. */

void
modularity_greedy_agglomerative(graph_t *g, char *alg_type,
				attr_id_t *membership,
				attr_id_t *num_communities,
				double *modularity)
{
    aggc_comm_t *communities; /* Size n to start off. The comm. adjacency lists
				 can be resized when necessary. */
    aggc_maxheap_t *maxheap; /* Max heap which stores community pairs that may
				result in the max modularity change if merged */

    aggc_adjcomm_t *comm_memchunk;
    double mod_val;
    long num_uniq_edges;
    aggc_heap_t *heap;
    attr_id_t curr_pair;
    attr_id_t *index;
    const attr_id_t *numEdges;
    const attr_id_t *endV;
    attr_id_t n, m, n_before, n_nz;
    attr_id_t i, j, v;
    attr_id_t start_iter, end_iter, new_start_iter;
    double max_key, new_max_key, m_inv, new_wt;
    attr_id_t max_key_idx, uniq_adj_count, v_degree_i, degree_i;
    aggc_adjcomm_t *comm_adjcomm, *adj_buffer;
    attr_id_t no_of_joins, total_joins;
    attr_id_t final_num_communities;
    double rij;
    double total_weight;
    weight_t L;
#if PREPROC_BICONN
    attr_id_t *bcc_num;
    attr_id_t *bcc_counts;
    attr_id_t bcc_total;
#endif

    printf("alg_type: %s\n",alg_type);

    keytype = NUMKEYTYPE;

    if (!strcmp(alg_type,"CNM"))
	keytype = KEY_CNM;
    else
	if (!strcmp(alg_type,"MB"))
	    keytype = KEY_MB;
	else
	    if (!strcmp(alg_type,"RAT"))
		keytype = KEY_RAT;
	    else
		if (!strcmp(alg_type,"MBRAT"))
		    keytype = KEY_MBRAT;
		else
		    if (!strcmp(alg_type,"LIN"))
			keytype = KEY_LIN;
		    else
			fprintf(stderr,"ERROR: unknown Algorithm: %s\n",alg_type);

    assert(keytype < NUMKEYTYPE);


#if PREPROC_BICONN
    bcc_num = (attr_id_t *)malloc(g->n*sizeof(attr_id_t));
    assert(bcc_num != NULL);
    bcc_total = biconnected_components(g, bcc_num);
    fprintf(stdout,"Biconnected Components: %12d\n",bcc_total);

    bcc_counts = (attr_id_t *)malloc(g->n*sizeof(attr_id_t));
    assert(bcc_counts != NULL);
    for (i=0 ; i<g->n ; i++)
	bcc_counts[i] = 0;
#endif

    n = g->n;
    m = g->m;        /* we store each undirected edge twice,
			and so m = 2 * the actual number of edges */
    numEdges = g->numEdges;
    endV = g->endV;
    mod_val = 0.0;

    switch (g->weight_type) {
	case 0: /* unweighted */
	    L.weighted = 0;
	    break;
	case 1: /* integer weights */
	case 2: /* long weights */
	case 3: /* float weights */
	case 4: /* double weights */
	    L.weighted = 1;
	    break;
	default:
	    fprintf(stderr,"ERROR: This weight type (%d) for modularity_greedy_agglomerative not implemented\n",
		    g->weight_type);
	    exit(-1);
    }

    /* Initialize */

    comm_memchunk = (aggc_adjcomm_t *)malloc(m * sizeof(aggc_adjcomm_t));
    assert(comm_memchunk != NULL);
    adj_buffer = (aggc_adjcomm_t *)malloc(n * sizeof(aggc_adjcomm_t));
    assert(adj_buffer != NULL);
    communities = (aggc_comm_t *)malloc(n * sizeof(aggc_comm_t));
    assert(communities != NULL);

    curr_pair = 0;

    num_uniq_edges = 0;
    total_weight = 0;

    for (i=0; i<n; i++) {
	communities[i].a = 0;
	communities[i].parent_id = -1;
	communities[i].resized = 0;
	communities[i].comm_size = 1;
#ifdef DO_COMM_WT
	communities[i].comm_wt = 0;
#endif
	start_iter = numEdges[i];
	end_iter   = numEdges[i+1];
	v_degree_i   = end_iter - start_iter;
	communities[i].max_size = v_degree_i;
	for (j=start_iter; j<end_iter; j++) {
	    v = endV[j];
	    comm_memchunk[j].comm_id = v;
	    if (L.weighted) {

		switch (g->weight_type) {
		    case 1: /* integer weights */
			new_wt = (double)g->int_weight_e[j];
			break;
		    case 2: /* long weights */
			new_wt = (double)g->l_weight_e[j];
			break;
		    case 3: /* float weights */
			new_wt = (double)g->fl_weight_e[j];
			break;
		    case 4: /* double weights */
			new_wt = g->dbl_weight_e[j];
			break;
		    default:
			fprintf(stderr,"ERROR: This weight type (%d) for modularity_greedy_agglomerative not implemented\n",
				g->weight_type);
			exit(-1);
		}
		comm_memchunk[j].Xij  = new_wt;
		total_weight         += new_wt;
		communities[i].a     += new_wt;
	    }
	    else {
		comm_memchunk[j].lij  = 1;
	    }
	}
	qsort(comm_memchunk+start_iter, v_degree_i, sizeof(aggc_adjcomm_t), comp_comm);
	uniq_adj_count = 0;
	if (v_degree_i > 0) {
	    if (comm_memchunk[start_iter].comm_id != i) {
		uniq_adj_count = 1;
		for (j=start_iter+1; j<end_iter; j++) {
		    if (comm_memchunk[j].comm_id != i) {
			if (comm_memchunk[j].comm_id !=
				comm_memchunk[start_iter+uniq_adj_count-1].comm_id) {
#if 0
			    comm_memchunk[start_iter+uniq_adj_count++].comm_id
				= comm_memchunk[j].comm_id;
#else
			    memcpy(comm_memchunk + start_iter + uniq_adj_count, comm_memchunk + j, sizeof(aggc_adjcomm_t));
			    uniq_adj_count++;
#endif
			}
		    }
		}
	    } else {
		uniq_adj_count = 0;
		for (j=start_iter+1; j<end_iter; j++) {
		    if (comm_memchunk[j].comm_id != i) {
#if 0
			comm_memchunk[start_iter+uniq_adj_count++].comm_id
			    = comm_memchunk[j].comm_id;
#else
			memcpy(comm_memchunk + start_iter + uniq_adj_count, comm_memchunk + j, sizeof(aggc_adjcomm_t));
			uniq_adj_count++;
#endif
			new_start_iter = j;
			break;
		    }
		}
		if (j < end_iter) {
		    uniq_adj_count = 1;
		    for (j=new_start_iter+1; j<end_iter; j++) {
			if (comm_memchunk[j].comm_id != i) {
			    if (comm_memchunk[j].comm_id !=
				    comm_memchunk[start_iter+uniq_adj_count-1].comm_id) {

#if 0
				comm_memchunk[start_iter+uniq_adj_count++].comm_id
				    = comm_memchunk[j].comm_id;
#else
				memcpy(comm_memchunk + start_iter + uniq_adj_count, comm_memchunk + j, sizeof(aggc_adjcomm_t));
				uniq_adj_count++;
#endif
			    }
			}
		    }
		}
	    }
	}
	communities[i].adjcomm = comm_memchunk+start_iter;
	communities[i].degree = uniq_adj_count;
	if (!L.weighted) {
#if 0
	    fprintf(stderr,"communities[%8d].a: %9f  uniq_adj_count: %8d\n",i,communities[i].a, uniq_adj_count);
	    assert (communities[i].a == (double)uniq_adj_count);
#else
	    communities[i].a = (double)uniq_adj_count;
#endif
	}
	num_uniq_edges += uniq_adj_count;
    }


#if 1
    if ((keytype==KEY_LIN) && (!L.weighted)) {
	/* Make this a weighted graph! */
	L.weighted = 1;
	total_weight = (double)num_uniq_edges;

	for (i=0; i<n; i++) {
	    degree_i = communities[i].degree;
	    comm_adjcomm = communities[i].adjcomm;
	    for (j=0; j<degree_i; j++)
		comm_adjcomm[j].Xij = 1.0;
	}
    }
#endif


    if (L.weighted) {
	L.M = total_weight;
	L.origM = L.M;
#ifdef VERBOSE
	fprintf(stderr,"GRAPH WEIGHTED n: %d m: %d L.M: %f\n",n, m, L.M);
#endif
	m_inv  = 1.0/L.M;
    } else {
	L.A = num_uniq_edges;
#ifdef VERBOSE
	fprintf(stderr,"GRAPH_UNWEIGHTED n: %d m: %d L.A: %d\n",n, m, L.A);
#endif
	m_inv  = 1.0/(double)L.A;
    }

    for (i=0; i<n; i++)
	communities[i].a = communities[i].a * m_inv;

    if (keytype==KEY_LIN) {
	n_nz = 0;
	for (i=0 ; i<n ; i++) {
#ifdef DO_ALIN
	    communities[i].a_lin = communities[i].a;
#endif
	    if (communities[i].a > 0)
		n_nz++;
	}
    }

    for (i=0; i<n; i++) {
	degree_i = communities[i].degree;
	comm_adjcomm = communities[i].adjcomm;

	for (j=0; j<degree_i; j++) {
	    v = comm_adjcomm[j].comm_id;
	    comm_adjcomm[j].value[KEY_CNM] =
#if 0
		2.0 * (m_inv - (((double)degree_i * m_inv) * ((double)communities[v].degree * m_inv)) );
#else
	    (L.weighted ?
	     2.0 * ((comm_adjcomm[j].Xij * m_inv) - (communities[i].a * communities[v].a))
	     :
	     2.0 * (m_inv - (((double)degree_i * m_inv) * ((double)communities[v].degree * m_inv)))
	    );
#endif


#if 0
	    if (L.weighted)
		fprintf(stderr,"INITv %4d %4d  m_inv: %15.7f  V: %f n: %d Xij %15.7f c[i].a: %15.7f  c[v].a: %15.7f  vN: %15.7f v: %15.7f\n",
			i, v, m_inv, (2.0*total_weight)/((double)(n*n)), n,
			comm_adjcomm[j].Xij,
			communities[i].a,
			communities[v].a,
			((comm_adjcomm[j].Xij * m_inv) - (communities[i].a * communities[v].a)),
			comm_adjcomm[j].value[KEY_CNM]
		       );
#endif

	}

	max_key     = NO_MAX_KEY;
	max_key_idx = -1;
	for (j=0; j<degree_i; j++) {

	    if ((keytype == KEY_MB) || (keytype == KEY_RAT) || (keytype == KEY_MBRAT)) {
#if 0
		fprintf(stderr,"Init v1 L=num_uniq_edges/2 (%d)  qij=%15.10f  2L*qij: %7.6f sqrt(lij=1 - L*qij) (%9.6f)\n",
			num_uniq_edges/2,
			comm_adjcomm[j].value[KEY_CNM],
			(double)(num_uniq_edges) * comm_adjcomm[j].value[KEY_CNM],
			sqrt(1.0 - ((double)(num_uniq_edges/2) * comm_adjcomm[j].value[KEY_CNM])));
#endif
#if 0
		if (L.weighted)
		    fprintf(stderr,"Init v1 %4d %4d  L (%f)  qij=%15.10f  1/2*L*qij: %7.6f sqrt(lij=Xij - L*qij) (%9.6f)\n",
			    i, comm_adjcomm[j].comm_id,
			    L.M,
			    comm_adjcomm[j].value[KEY_CNM],
			    0.5 * L.M * comm_adjcomm[j].value[KEY_CNM],
			    sqrt(comm_adjcomm[j].Xij - (L.M * comm_adjcomm[j].value[KEY_CNM])));
#endif
		rij = 0.5 * (L.weighted ? L.M : (double)L.A) * comm_adjcomm[j].value[KEY_CNM];
#if 1
		if (L.weighted)
		    assert (rij < comm_adjcomm[j].Xij);
		else
		    assert (rij < 1.0);
#endif
		comm_adjcomm[j].value[KEY_MB] =  rij / sqrt( (L.weighted ? comm_adjcomm[j].Xij : 1.0) - rij);
	    }

	    if (keytype != KEY_CNM) {
		/* Initial cratio = 1, since each community has size 1 */
		comm_adjcomm[j].value[KEY_RAT]   =  comm_adjcomm[j].value[KEY_CNM];
		comm_adjcomm[j].value[KEY_MBRAT] =  comm_adjcomm[j].value[KEY_MB];

		if (keytype == KEY_LIN) {

		    assert(L.weighted);

		    v = comm_adjcomm[j].comm_id;

		    /* val = Xij - (V + V1(i) + V2(j)), where V = mean, and V1(i) = sum/wt - V */
		    comm_adjcomm[j].value[KEY_LIN] =
			2.0 *
			(comm_adjcomm[j].Xij - (L.M * (communities[i].a + communities[v].a)/(double)n_nz) + (L.M /(double)(n_nz*n_nz)) );
#ifdef DEBUG_LIN
		    fprintf(stderr,"DEBUG_LIN Init v1 %4d %4d  L (%f) Xij: %15f i.a: %8f j.a: %8f n: %d n2: %d v: %15f \n",
			    i, v,
			    L.M,
			    comm_adjcomm[j].Xij,
			    communities[i].a,
			    communities[v].a,
			    n_nz,
			    n_nz*n_nz,
			    comm_adjcomm[j].value[KEY_LIN]);
#endif
		}

#ifdef CULL_DEGREE_ONE
		if ((degree_i == 1) || (communities[comm_adjcomm[j].comm_id].degree == 1))
		    comm_adjcomm[j].value[keytype] = INFTY;
#endif
#if PREPROC_BICONN
#if 1
		fprintf(stderr,"BICONN edge %12d (%4d), %12d (%4d)\n",
			i, bcc_num[i], comm_adjcomm[j].comm_id, bcc_num[comm_adjcomm[j].comm_id]);
		if (bcc_num[i] != bcc_num[comm_adjcomm[j].comm_id]) {
		    fprintf(stderr,"BICONN BRIDGE %12d (%4d), %12d (%4d)\n",
			    i, bcc_num[i], comm_adjcomm[j].comm_id, bcc_num[comm_adjcomm[j].comm_id]);
		    bcc_counts[i]++;
		    bcc_counts[comm_adjcomm[j].comm_id]++;
		}
#endif
		if (bcc_num[i] != bcc_num[comm_adjcomm[j].comm_id])
		    comm_adjcomm[j].value[keytype] = 0;
#endif
	    }

	    if (comm_adjcomm[j].value[keytype] > max_key) {
		max_key     = comm_adjcomm[j].value[keytype];
		max_key_idx = j;
	    }
	}
	communities[i].max_key     = max_key;
	communities[i].max_key_idx = max_key_idx;
#if PREPROC_BICONN
	fprintf(stderr,"BICONN vert %12d %4d\n", i, bcc_num[i]);
#endif
    }

#if PREPROC_BICONN
    for (i=0 ; i<n ; i++)
	if (bcc_counts[i] > 4)
	    fprintf(stdout,"BICONN T NODE: %12d  (%12d)\n",i, bcc_counts[i]);
#endif

    /* DEBUG */

    fprintf(stderr, "Number of edges after removing self loops and duplicates: %ld\n",
	    num_uniq_edges/2);

#ifdef DEBUG
    aggc_print_communities(communities, n, L);
#endif

    /* Compute initial modularity */
    mod_val = 0.0;
    for (i=0; i<n; i++) {
	mod_val -= communities[i].a * communities[i].a;
    }

    maxheap = (aggc_maxheap_t *) malloc(sizeof(aggc_maxheap_t));
    assert(maxheap != NULL);

    maxheap->heap = (aggc_heap_t *) malloc((n+1) * sizeof(aggc_heap_t));
    assert(maxheap->heap != NULL);
    maxheap->index = (attr_id_t *) malloc(n * sizeof(attr_id_t));
    assert(maxheap->index != NULL);

    heap  = maxheap->heap;
    index = maxheap->index;

    maxheap->n = 0;

#ifdef DO_HEAP2
    if (keytype==KEY_LIN)
	maxheap->use2 = 1;
    else
	maxheap->use2 = 0;
#endif

    /* heapify */
    for (i=0; i<n; i++) {
	if (communities[i].max_key == NO_MAX_KEY) {
	    index[i] = -1;
	    continue;
	}
	index[i]    = maxheap->n;
	max_key     = communities[i].max_key;
	max_key_idx = communities[i].max_key_idx;
	heap[maxheap->n].comm_id = i;
	heap[maxheap->n].val     = communities[i].adjcomm[max_key_idx].value[keytype];
#ifdef DO_HEAP2
	if (maxheap->use2)
	    heap[maxheap->n].val2  = communities[i].adjcomm[max_key_idx].value[0];
#endif
	aggc_maxheap_sift_up(maxheap, maxheap->n);
	maxheap->n++;
    }

#if 1
    fprintf(stderr, "n: %d   Number of initial communities: %d\n", n, maxheap->n);
#endif

    aggc_maxheap_check(maxheap);

    for (MaxKeyHistoryIdx = 0 ; MaxKeyHistoryIdx < MAXHISTORY ; MaxKeyHistoryIdx++)
	MaxKeyHistory[MaxKeyHistoryIdx] = 0.0;
    MaxKeyHistoryIdx    = 0;
    MaxKeyHistoryFilled = 0;
    sumMKH  = 0.0;
    sumMKH2 = 0.0;

#ifdef DEBUG_MAXHEAP
    for (i=0; i<maxheap->n; i++) {
	fprintf(stderr, "DEBUG_MAXHEAP maxheap->heap[%4d].val (%20f) "
#ifdef DO_HEAP2
		"val2 (%20f) "
#endif
		"comm_id (%4d)\n",
		i, maxheap->heap[i].val,
#ifdef DO_HEAP2
		maxheap->heap[i].val2,
#endif
		maxheap->heap[i].comm_id);
    }

    for (i=0; i<maxheap->n; i++) {
	fprintf(stderr, "maxheap->heap[index[%4d] (%4d)].comm_id: %4d\n", i, index[i], maxheap->heap[index[i]].comm_id);
    }
#endif

    total_joins = maxheap->n - 1;

    no_of_joins = 0;

#ifdef VERBOSE
    fprintf(stderr, "Initial modularity: %lf\n", mod_val);
#endif

    while (no_of_joins < total_joins) {

	n_before = maxheap->n;

#if 0
	for (i=0; i<maxheap->n; i++)
	    for (j=0 ; j<communities[i].degree ; j++)
		fprintf(stderr, "bef merge %4d %4d %9.6f %9.6f %9.6f\n",
			i,
			communities[i].adjcomm[j].comm_id,
			communities[i].adjcomm[j].value[KEY_CNM],
			communities[i].adjcomm[j].value[KEY_MB],
			communities[i].adjcomm[j].value[KEY_CNM] / communities[i].adjcomm[j].value[KEY_MB]);
#endif

#ifdef DEBUG_MAXKEY
	aggc_maxkey_check(communities,maxheap);
#endif

	new_max_key = aggc_merge_communities(communities, maxheap, adj_buffer, n - no_of_joins, &L);

#if 0
	final_num_communities = aggc_compute_membership(communities, membership, n);
	for (i=0 ; i<n ; i++)
	    if (g->zero_indexed)
		fprintf(stdout,"COMMUNITY %12d %12d %12d\n",no_of_joins, i, membership[i]+1);
	    else
		fprintf(stdout,"COMMUNITY %12d %12d %12d\n",no_of_joins, i+1, membership[i]+1);
#endif

#ifdef DEBUG_MAXHEAP
	for (i=0; i<maxheap->n; i++) {
	    fprintf(stderr, "DEBUG_MAXHEAP maxheap->heap[%4d].val (%20f) "
#ifdef DO_HEAP2
		    "val2 (%20f) "
#endif
		    "comm_id (%4d)\n",
		    i, maxheap->heap[i].val,
#ifdef DO_HEAP2
		    maxheap->heap[i].val2,
#endif
		    maxheap->heap[i].comm_id);
	}
#endif


#ifdef DEBUG_MAXHEAP
	aggc_maxheap_check(maxheap);
#endif

#ifdef DEBUG
	aggc_print_communities(communities, n, L);
#endif

	no_of_joins++;

#if 0
	if (new_max_key < 0) {
	    break;
	}
#else
	if (maxheap->n == n_before) {
	    break;
	}
#endif

	mod_val += new_max_key;
	if ((no_of_joins % 10000) == 0)
	    fprintf(stderr, "join: %d, mod %9.6f\n", no_of_joins, mod_val);
    }

#if 1
    fprintf(stderr, "n: %d   Number of non-singleton communities: %d\n",
	    n, maxheap->n);
#endif

#if 0
    fprintf(stderr, "\n");
    /*  aggc_print_communities(communities, n, L); */
    aggc_print_final_communities(communities, maxheap, L);
#endif

#if 1
    final_num_communities = aggc_compute_membership(communities, membership, n);
#else

    /* get final community membership information */
    final_num_communities = 0;
    for (i=0; i<n; i++)
	if (communities[i].parent_id == -1)
	    membership[i] = final_num_communities++;

    for (i=0; i<n; i++) {
	if (communities[i].parent_id != -1) {
	    parent_id = communities[i].parent_id;
	    while (communities[parent_id].parent_id != -1)
		parent_id = communities[parent_id].parent_id;
	    membership[i] = membership[parent_id];
	}
    }
#endif

    *num_communities = final_num_communities;
    *modularity = mod_val;

    for (i=0; i<n; i++)
	if (communities[i].resized == 1)
	    free(communities[i].adjcomm);

    free(adj_buffer);
    free(comm_memchunk);
    free(communities);
    free(maxheap->heap);
    free(maxheap->index);
    free(maxheap);
#if PREPROC_BICONN
    free(bcc_counts);
    free(bcc_num);
#endif

    return;
}
