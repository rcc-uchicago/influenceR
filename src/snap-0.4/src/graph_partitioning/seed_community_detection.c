#include <string.h>
#include <time.h>
#include "graph_partitioning.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include <math.h>

#include "internal-timer.h"
#include "modularity_greedy_agglomerative.h"

/** Greedily grow a region around the provided seeds until the
    modularity would decrease.  The style of growth is controlled by
    alg_type.

 Possible algorithms:
   - "CNM" : CNM-style local maximization
   - "MB" : McCloskey-Bader normalization and relevance tests
   - "LIN" : Linear background model
   - "RAT", "MBRAT" : Ratio tests for CNM and MB.

 @param g Input graph (routine only tested for undirected)
 @param alg_type Algorithm to use
 @param seeds Array of seed vertex indices
 @param n_seeds Length of seeds
 @param membership Output membership array, -1 for vertices not in the community.
 @param modularity Accumulated modularity for grown region
*/
void
seed_set_community_detection(graph_t *g, char *alg_type,
			     attr_id_t* seeds, attr_id_t n_seeds,
			     attr_id_t *membership,
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

  /*printf("alg_type: %s\n",alg_type);*/

  tic ();

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
	      memmove(comm_memchunk + start_iter + uniq_adj_count, comm_memchunk + j, sizeof(aggc_adjcomm_t));
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
	    memmove(comm_memchunk + start_iter + uniq_adj_count, comm_memchunk + j, sizeof(aggc_adjcomm_t));
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
		memmove(comm_memchunk + start_iter + uniq_adj_count, comm_memchunk + j, sizeof(aggc_adjcomm_t));
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

#ifdef DEBUG
  fprintf(stderr, "Number of edges after removing self loops and duplicates: %ld\n", 
	  num_uniq_edges/2);

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
  for (attr_id_t k = 0; k < n; ++k)
    index[k] = -1;
  for (attr_id_t seedk = 0; seedk < n_seeds; ++seedk) {
    i = seeds[seedk];
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

#if 0
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

  total_joins = n - 1;

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

#if 1
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

#if 0
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

  toc ();

  return;
}

static inline int
approx_pagerank_push(const graph_t *g, attr_id_t u, double * p, double *r, double alpha,
		     int * queue, int end, double epsilon)
{	
  attr_id_t start = g->numEdges[u];
  attr_id_t end2 = g->numEdges[u+1];
  attr_id_t deg = end2 - start;

  int count = 0;

  attr_id_t i;
  //printf("starting loop %d to %d\n", start, end);
  while(r[u]/(double) deg >= epsilon)
    {
      for(i=start; i<end2; i++)
	{
	  attr_id_t v = g->endV[i];
	  //fprintf(stderr, "changing i=%d from %lf\n", i, r[i]);
	  double rprev = r[v];
	  r[v] = r[v] + (1-alpha)*r[u]/(2.0*deg);

	  if(rprev < epsilon && r[v] >= epsilon)
	    {
	      queue[end+count] = v;
	      count++;
	    }
	}
      //fprintf(stderr, "old ru = %lf, pu = %lf, u = %d\n", r[u], p[u], u);
      r[u] = (1-alpha)*r[u]/2.0;
      p[u] += alpha*r[u];
    }
  //fprintf(stderr, "new ru = %lf, pu = %lf, u = %d\n", r[u], p[u], u);

  return count;
}

struct vp {
  attr_id_t v;
  attr_id_t deg;
  double p;
  double weighted_p;
};

static int
compare_vp (const void * a, const void * b)
{
  double r = ((struct vp*)a)->weighted_p - ((struct vp*)b)->weighted_p;
  if (r > 0) return -1;
  if (r < 0) return 1;
  return 0; /* NaNs fall to here. */
}

static void
compute_approx_pagerank (const graph_t *g, const attr_id_t * seeds, int nseed,
			 double alpha, double eps,
			 double *p)
{
  const attr_id_t N = g->n;
  const attr_id_t * restrict xoff = g->numEdges;

  double * r = malloc(sizeof(double)*N + sizeof(int)*g->m);
  int * queue = (int *) &r[N];

  int start, end;
  attr_id_t i;

  for(i=0; i<N; i++)
    {
      p[i] = 0;
      r[i] = 0;
    }
  for(i=0; i<nseed; i++)
    r[seeds[i]] = 1.0/nseed;

  for(i=0; i<nseed; i++)
    queue[i] = seeds[i];
  start = 0;
  end = nseed;

  while(start < end)
    {
      attr_id_t deg;
      i = queue[start];
      deg = xoff[i+1] - xoff[i];

      if(deg != 0 && r[i]/(double) deg >= eps)
	end += approx_pagerank_push (g, i, p, r, alpha, queue, end, eps);
      start++;
    }

  free (r);
}

/** Grow a region around the provided seeds using a personalized
    PageRank algorithm and conductance minimization.

 @param g Input graph (routine only tested for undirected)
 @param seeds Array of seed vertex indices
 @param n_seeds Length of seeds
 @param alpha Regularization parameter for PageRank
 @param eps Approximation factor for approximate PageRank
 @param membership_in Output membership array, -1 for vertices not in the community.
*/
void
pagerank_community (graph_t *g, const attr_id_t * seed, attr_id_t nseed,
		    attr_id_t * membership_in,
		    double alpha, double eps)
{
  attr_id_t * restrict membership = membership_in;
  const attr_id_t N = g->n;
  const attr_id_t * restrict xoff = g->numEdges;
  const attr_id_t * restrict xadj = g->endV;
  const int nseed_seen = 0;

  tic ();

  struct vp *vlist = malloc (N * (sizeof (*vlist) + sizeof (double) + sizeof (int)));
  double *p = (double*)&vlist[N];
  int *loc = (int*)&p[N];

  double mincond = HUGE_VAL;
  attr_id_t mink = -1, vol_s, L_s;

  compute_approx_pagerank (g, seed, nseed, alpha, eps, p);

  for (attr_id_t k = 0; k < N; ++k) {
    vlist[k].v = k;
    vlist[k].weighted_p = p[k];
  }
  qsort (vlist, N, sizeof(*vlist), compare_vp);
  for (attr_id_t k = 0; k < N; ++k)
    loc[vlist[k].v] = k;

  /* Compute the running conductance, saving the least. */
  vol_s = L_s = 0;
  for (attr_id_t k = 0; k < N; ++k) {
    const attr_id_t v = vlist[k].v;
    const attr_id_t deg =  xoff[v+1] - xoff[v];
    double cond;

    vol_s += deg;
    for (attr_id_t k2 = xoff[v]; k2 < xoff[v+1]; ++k2) {
      const attr_id_t w = xadj[k2];
      if (w == v) continue;
      if (loc[w] >= 0 && loc[w] < k)
	++L_s; /* Will only count internal edges *once*. */
    }
    {
      double d_s, vol_bar_s, minvol;
      d_s = vol_s - 2*L_s;
      vol_bar_s = g->m - 2*L_s;
      minvol = (vol_s < vol_bar_s? vol_s : vol_bar_s);
      if (minvol)
	cond = d_s / minvol;
      else
	cond = HUGE_VAL;
      /* fprintf (stderr, "v %ld vol_s %ld L_s %ld d_s %ld vol_bar_s %ld cond %g\n", */
      /* 	       (long)v, */
      /* 	       (long)vol_s, (long)L_s, (long)d_s, (long)vol_bar_s, */
      /* 	       cond); */
      assert (g->m == vol_s + vol_bar_s - d_s);
    }
    assert (vol_s >= 2*L_s);
    assert (cond >= 0.0);

    if (cond < mincond) {
      mincond = cond;
      mink = k;
    }
  }

  /* Ensure all seeds are included. */
  for (attr_id_t k = 0; k < nseed; ++k)
    if (mink < loc[seed[k]]) mink = loc[seed[k]];

  /* Save the global best to the output array, membership[]. */
  memset (membership, 0, N * sizeof (*membership));
  for (attr_id_t k = 0; k < mink; ++k) {
    const attr_id_t v = vlist[k].v;
    /* Do not include degree-0 vertices. */
    if (0 != xoff[v+1] - xoff[v])
      membership[v] = 1;
  }

  toc ();
}

/** Grow a region around the provided seeds using a limited version
    of Andersen and Lang's random walk algorithm.

 @param g Input graph (routine only tested for undirected)
 @param seeds Array of seed vertex indices
 @param n_seeds Length of seeds
 @param membership_in Output membership array, -1 for vertices not in the community.
 @param nsweep Number of sweeps
*/
void
andersen_lang (graph_t *g, const attr_id_t * seed, attr_id_t nseed,
	       attr_id_t * membership_in, const int nsweep)
{
  attr_id_t * restrict membership = membership_in;
  const attr_id_t N = g->n;
  const attr_id_t * restrict xoff = g->numEdges;
  const attr_id_t * restrict xadj = g->endV;

  tic ();

  attr_id_t nv = 0;
  struct vp *vlist = malloc (N * (sizeof (*vlist) + sizeof (attr_id_t) + sizeof (int)));
  attr_id_t *best_so_far = (attr_id_t*)&vlist[N];
  int *loc = (int*)&best_so_far[N];
  double sum;
  int sweepno;
  double global_mincond = HUGE_VAL;
  attr_id_t saved_nv;

  for (attr_id_t k = 0; k < N; ++k)
    loc[k] = -1;

  /* Initialize with the seeds. */
  sum = 0.0;
  for (attr_id_t k = 0; k < nseed; ++k) {
    const attr_id_t v = seed[k];
    const attr_id_t deg = xoff[v+1] - xoff[v];
    vlist[nv].v = v;
    vlist[nv].deg = deg;
    vlist[nv].p = deg;
    sum += deg;
    loc[v] = nv;
    ++nv;
  }
  for (attr_id_t k = 0; k < nseed; ++k) {
    vlist[nv].p /= sum;
    vlist[nv].weighted_p = vlist[nv].p;
  }


  for (int sweepno = 0; sweepno < nsweep; ++sweepno) {
    const attr_id_t NV = nv;
    attr_id_t vol_s, L_s;
    double mincond = HUGE_VAL;
    attr_id_t mink = 0, minseedk = -1;
    /* In-place sparse matrix - sparse vector product: vlist = 1/2 (I + D\A) * vlist */
    for (attr_id_t k = 0; k < NV; ++k) {
      loc[vlist[k].v] = k; /* Can be over-written each time: no deletions. */
      vlist[k].weighted_p = vlist[k].p; /* temporarily hold old value */
      vlist[k].p = 0.5 * vlist[k].p; /* Accumulate 1/2 * I . */
    }
    for (attr_id_t k = 0; k < NV; ++k) {
      const attr_id_t v = vlist[k].v;
      for (attr_id_t k2 = xoff[k]; k2 < xoff[k+1]; ++k2) {
	const attr_id_t w = xadj[k2];
	const attr_id_t degw = xoff[w+1] - xoff[w];
	const double x = vlist[k].weighted_p; /* Old value. */
	assert (degw > 0);
	if (w == v) continue; /* Skip self-edges for no terribly good reason. */
	if (loc[w] < 0) { /* New vertex. */
	  loc[w] = ++nv;
	  vlist[loc[w]].v = w;
	  vlist[loc[w]].deg = degw;
	  vlist[loc[w]].p = 0.0;
	}
	vlist[loc[w]].p += 0.5 * x / degw; /* Accumulate 1/2 D\A . */
      }
    }
    for (attr_id_t k = 0; k < nv; ++k)
      vlist[k].weighted_p = (vlist[k].deg? vlist[k].p / vlist[k].deg : 0.0);

    /* Sort, relocate. */
    qsort (vlist, nv, sizeof(*vlist), compare_vp);
    for (attr_id_t k = 0; k < nv; ++k)
      loc[vlist[k].v] = k;

    /* Ensure the seeds are included. */
    for (attr_id_t k = 0; k < nseed; ++k) {
      assert(loc[seed[k]] >= 0);
      if (loc[seed[k]] > minseedk)
	minseedk = loc[seed[k]];
    }

    /* Compute the running conductance, saving the least. */
    vol_s = L_s = 0;
    for (attr_id_t k = 0; k < nv; ++k) {
      const attr_id_t v = vlist[k].v;
      const attr_id_t deg = xoff[v+1] - xoff[v];
      double cond;

      vol_s += deg;
      for (attr_id_t k2 = xoff[v]; k2 < xoff[v+1]; ++k2) {
	const attr_id_t w = xadj[k2];
	if (w == v) continue;
	if (loc[w] >= 0 && loc[w] < k)
	  ++L_s; /* Will only count internal edges *once*. */
      }
      {
	double d_s, vol_bar_s, minv;
	d_s = vol_s - 2*L_s;
	vol_bar_s = g->m - 2*L_s;
	minv = (vol_s < vol_bar_s? vol_s : vol_bar_s);
	if (minv > 0)
	  cond = d_s / minv;
	else if (d_s == 0)
	  cond = 1.0;
	else
	  cond = HUGE_VAL;
	/* fprintf (stderr, "%ld: vol_s %ld L_s %ld d_s %ld vol_bar_s %ld cond %g\n", */
	/* 	 (long)sweepno, */
	/* 	 (long)vol_s, (long)L_s, (long)d_s, (long)vol_bar_s, */
	/* 	 cond); */
	assert (g->m == vol_s + vol_bar_s - d_s);
      }
      assert (vol_s >= 2*L_s);
      assert (cond >= 0.0);

      /* All seeds must be included in every iteration. */
      if (cond < mincond && k >= minseedk) {
	mincond = cond;
	mink = k;
      }
    }

    /* Copy the least to the temporary output if less than the global least. */
    if (mincond < global_mincond) {
      global_mincond = mincond;
      saved_nv = mink;
      for (attr_id_t k = 0; k <= mink; ++k) {
	best_so_far[k] = vlist[k].v;
      }
    }
  }

  /* Save the global best to the output array, membership[]. */
  memset (membership, 0, N * sizeof (*membership));
  for (attr_id_t k = 0; k < saved_nv; ++k) {
    const attr_id_t v = best_so_far[k];
    /* Do not include degree-0 vertices. */
    if (0 != xoff[v+1] - xoff[v])
      membership[v] = 1;
  }

  toc ();
}


/** Grow a region around the provided seeds with a fixed number of
    breadth-first steps.

 @param g Input graph (routine only tested for undirected)
 @param seeds Array of seed vertex indices
 @param num_seeds Length of seeds
 @param membership_in Output membership array, -1 for vertices not in the community.
 @param steps Number of steps
*/
void
BFS_seed_set_expansion (graph_t *g, attr_id_t *seeds, int num_seeds, attr_id_t *membership,
			attr_id_t steps)
{
  const attr_id_t N = g->n;
  const attr_id_t * restrict xoff = g->numEdges;
  const attr_id_t * restrict xadj = g->endV;
  attr_id_t * restrict queue = NULL;
  attr_id_t qi, qj, step;

  tic ();

  memset (membership, 0, N * sizeof (*membership));

  queue = malloc (N * sizeof (*queue));
  qi = 0;
  for (qj = 0; qj < num_seeds; ++qj) {
    const attr_id_t v = seeds[qj];
    queue[qj] = v;
    membership[v] = 1;
  }

  for (attr_id_t step = 0; step < steps; ++step) {
    const attr_id_t qend = qj;
    /*fprintf (stderr, "BFS step %d, len %ld\n", (int)step, (long)(qend-qi));*/
    for (; qi < qend; ++qi) {
      const attr_id_t v = queue[qi];
      for (attr_id_t k = xoff[v]; k < xoff[v+1]; ++k) {
	const attr_id_t w = xadj[k];
	if (!membership[w]) {
	  membership[w] = 1;
	  queue[qj++] = w;
	}
      }
    }
    if (qend == qj) break;
  }

  free (queue);
  toc ();
}

static attr_id_t
random_walk(graph_t *g, attr_id_t V, int num_levels)
{
  int i;
  attr_id_t index, deg_V, rand_num;

  for(i=0; i<num_levels; i++)
    {

      // number of vertices adjacent to the vertex V i.e. degree[V]
      deg_V = g->numEdges[V + 1] - g->numEdges[V];  

      if(deg_V == 0)                 //if degree of V is zero then new random vertex will be taken as the source
	{
	  return -2;
	}
      rand_num = (attr_id_t) (drand48() * deg_V);  // random number between 1 and deg_V

      index = g->numEdges[V] + rand_num;   //this will refer to an edge whose end vertex will be a random neighbor of vertex V

      V = g->endV[index];    //random neighbor of random_vertex
    }

  return(V);

}
