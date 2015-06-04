#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "graph_partitioning.h"
#include "utils.h"
#include "stdlib.h"

#include <float.h>
#include <time.h>

#if defined(_OPENMP)
#define OMP(x) _Pragma(x)
#if defined(__GNUC__)||defined(__INTEL_COMPILER)
/* The following assume the Intel / gcc spellings. */
#define atomic_val_compare_and_swap(out, a, v1, v2) \
  do { (out) = __sync_val_compare_and_swap ((a), (v1), (v2)); } while (0)
#define atomic_fetch_and_add(a, incr) __sync_fetch_and_add ((a), (incr))
#else
#error "Don't know spellings for atomic ops on this platform."
#endif
#elif defined(__MTA__)
/* Slightly bogus... Assumes integers. */
#define atomic_val_compare_and_swap(out, a, v1, v2)			\
  do { int xx_ = readfe ((a)); if (xx_ == v1) xx_ = v2;			\
    writeef ((a), xx_); (out) = xx_; } while (0)
#define atomic_fetch_and_add(a, incr) int_fetch_add ((a), (incr))
#else
#define OMP(x)
#define atomic_val_compare_and_swap(out, a, v1, v2)	\
  do { if (*(a) == (v1)) *(a) = (v2); (out) = *(a); } while (0)
#define atomic_fetch_and_add(a, incr) do { *(a) += (incr); } while (0)
#endif

#define INFTY (1<<30)

#define MAX_NUM_SEEDS 10

extern double snap_runtime;

static char *infilename;
static char *alg_type = NULL;
static graph_t *g = NULL;
static int num_communities;
static attr_id_t *membership = NULL, *member_flag = NULL;
static int round_robin = 1, stop_on_merge = 1, threshold = 0;
static double tol = 0.0;

static int num_seeds;
static attr_id_t seeds[MAX_NUM_SEEDS];
static attr_id_t comm_root[MAX_NUM_SEEDS];
static attr_id_t comm_size[MAX_NUM_SEEDS];
static attr_id_t comm_vol[MAX_NUM_SEEDS];
static attr_id_t comm_maxdeg[MAX_NUM_SEEDS];

static double modularity, clst, conductance;

static void
dump_output (int do_hdr)
{
  attr_id_t total_size = 0;
  attr_id_t total_vol = 0;
  attr_id_t total_maxdeg = 0;

  if (do_hdr) printf ("graphname");
  else printf ("%s", infilename);
  if (do_hdr) printf (", n");
  else printf (", %ld", (long)g->n);
  if (do_hdr) printf (", m");
  else printf (", %ld", (long)(g->m/2));

  if (do_hdr) printf (", alg_type");
  else printf (", %s", alg_type);
  if (do_hdr) printf (", RR");
  else printf (", %ld", (long)round_robin);
  if (do_hdr) printf (", stoponmerge");
  else printf (", %ld", (long)stop_on_merge);
  if (do_hdr) printf (", threshold");
  else printf (", %ld", (long)(threshold < INFTY? threshold : 0));
  if (do_hdr) printf (", tol");
  else printf (", %e", tol);

  if (do_hdr) printf (", nseed");
  else printf (", %ld", (long)num_seeds);
  for (int k = 0; k < MAX_NUM_SEEDS; ++k) {
    if (do_hdr) printf (", seed%ld", k);
    else printf (", %ld", (long)seeds[k]);
  }
  if (do_hdr) printf (", ncomm");
  else printf (", %ld", (long)num_communities);
  for (int k = 0; k < MAX_NUM_SEEDS; ++k) {
    if (do_hdr) printf (", comm_size%ld", k);
    else printf (", %ld", (long)comm_size[k]);
    total_size += (comm_size[k] >= 0? comm_size[k] : 0);
  }
  if (do_hdr) printf (", total_size");
  else printf (", %ld", (long)total_size);
  for (int k = 0; k < MAX_NUM_SEEDS; ++k) {
    if (do_hdr) printf (", comm_vol%ld", k);
    else printf (", %ld", (long)comm_vol[k]);
    total_vol += (comm_vol[k] >= 0? comm_vol[k] : 0);
  }
  if (do_hdr) printf (", total_vol");
  else printf (", %ld", (long)total_vol);
  for (int k = 0; k < MAX_NUM_SEEDS; ++k) {
    if (do_hdr) printf (", comm_maxdeg%ld", k);
    else printf (", %ld", (long)comm_maxdeg[k]);
    if (comm_maxdeg[k] > total_maxdeg) total_maxdeg = comm_maxdeg[k];
  }
  if (do_hdr) printf (", total_maxdeg");
  else printf (", %ld", (long)total_maxdeg);
  for (int k = 0; k < MAX_NUM_SEEDS; ++k) {
    if (do_hdr) printf (", comm_root%ld", k);
    else printf (", %ld", (long)comm_root[k]);
  }

  if (do_hdr) printf (", time");
  else printf (", %e", snap_runtime);
  if (do_hdr) printf (", modularity");
  else printf (", %e", modularity);
  if (do_hdr) printf (", conductance");
  else printf (", %e", conductance);
  if (do_hdr) printf (", clstcoeff");
  else printf (", %e", clst);

  printf ("\n");

  /*
  if (do_hdr) printf ("");
  else printf ("");
  */
}

static int compare_attr_id_t (const void * a, const void * b)
{
  return (int)( *(const attr_id_t*)a - *(const attr_id_t*)b );
}

static void
identify_comm (int *num_communities, attr_id_t *comm_root, attr_id_t *comm_size,
	       attr_id_t *membership, attr_id_t *member_lbl)
/*
  Given a vague labeling in member_lbl, identify communities by their
  first seed.  All vertices not in a community are labeled with -1.
*/
{
  const int N = g->n;
  int nc = 0;

  int ncomm;
  attr_id_t seed_label[MAX_NUM_SEEDS];
  graph_t *subG = NULL;
  attr_id_t *queue = NULL;
  attr_id_t *flag = NULL;

  const attr_id_t * restrict xoff = g->numEdges;
  const attr_id_t * restrict xadj = g->endV;

  queue = malloc (2 * N * sizeof (*queue));
  assert (queue);
  flag = &queue[N];

  memset (comm_size, 0, sizeof (comm_size));
  memset (comm_vol, 0, sizeof (comm_vol));
  memset (comm_maxdeg, 0, sizeof (comm_maxdeg));

  /* Find the seeds' labels, then convert member_lbl to an indicator array. */
  for (int k = 0; k < num_seeds; ++k) {
    seed_label[k] = member_lbl[seeds[k]];
    /* Ugh. CNM/MB aren't mergining neighboring communities. */
    /*fprintf (stderr, "seed label %ld = %ld\n", (long)k, (long)seed_label[k]);*/
  }

  OMP("omp parallel for")
  for (int i = 0; i < N; ++i) {
    const attr_id_t cid = member_lbl[i];
    membership[i] = -1;
    for (int k = 0; k < num_seeds; ++k)
      if (cid == seed_label[k])
	membership[i] = -2;
  }

  /* Find and label the connected components. */
  ncomm = 0;
  for (attr_id_t kseed = 0; kseed < num_seeds; ++kseed) {
    attr_id_t k1, k2;
    const attr_id_t seedv = seeds[kseed];
    if (membership[seedv] != -2) continue;
    queue[0] = seedv;
    membership[seedv] = seedv;
    comm_root[ncomm] = seedv;

    k1 = 0;
    k2 = 1;
    do {
      const attr_id_t qkend = k2;
      attr_id_t cv = 0;
      OMP("omp parallel for reduction(+:cv)")
      for (attr_id_t k1i = k1; k1i < qkend; ++k1i) {
	const attr_id_t v = queue[k1i];
	const attr_id_t deg = xoff[v+1] - xoff[v];
	cv += deg;
	if (deg > comm_maxdeg[ncomm]) comm_maxdeg[ncomm] = deg;
	for (attr_id_t k = xoff[v]; k < xoff[v+1]; ++k) {
	  const attr_id_t w = xadj[k];
	  attr_id_t memb;
	  if (membership[w] < -1) {
	    atomic_val_compare_and_swap (memb, &membership[w], -2, seedv);
	    if (memb < -1) {
	      //if (membership[w] < -1) {
	      //membership[w] = seedv;
	      //int loc;
	      //OMP("omp atomic") loc = k2++;
	      queue[atomic_fetch_and_add (&k2, 1)] = w;
	    }
	  }
	}
      }
      k1 = qkend;
      comm_vol[ncomm] = cv;
    } while (k1 != k2);
    comm_size[ncomm] = k2;
    ++ncomm;
  }
  for (attr_id_t k = 0; k < N; ++k)
    member_lbl[k] = (membership[k] < 0? 0 : 1);
  *num_communities = ncomm;
}

static void
run_MB_CNM (int *num_communities, attr_id_t *membership, attr_id_t *member_flag)
{
  double mod;
  attr_id_t inseeds[MAX_NUM_SEEDS];

  for (int i = 0; i < num_seeds; ++i) inseeds[i] = seeds[i];

  OMP("omp parallel for")
    for (int i = 0; i < g->n; ++i)
      member_flag[i] = -9999;

  seed_set_community_detection (g, alg_type, inseeds, num_seeds,
				member_flag, &mod);
  if (getenv ("VERBOSE")) {
    for (int i = 0; i < num_seeds; ++i)
      if (inseeds[i] != seeds[i])
	fprintf (stderr, "MBCNM %s altered seed %d, %d->%d\n",
		 alg_type, i, seeds[i], inseeds[i]);
    fprintf (stderr, "mod: %e\n", mod);
  }
  identify_comm (num_communities, comm_root, comm_size, membership, member_flag);
}

static void
run_BFS (int *num_communities, attr_id_t *membership, attr_id_t *member_flag)
{
  attr_id_t inseeds[MAX_NUM_SEEDS];

  for (int i = 0; i < num_seeds; ++i) inseeds[i] = seeds[i];

OMP("omp parallel for")
  for (int i = 0; i < g->n; ++i)
    member_flag[i] = -9999;

#define BFS_STEPS 3
  BFS_seed_set_expansion (g, inseeds, num_seeds, member_flag, BFS_STEPS);

  if (getenv ("VERBOSE"))
    for (int i = 0; i < num_seeds; ++i)
      if (inseeds[i] != seeds[i])
	fprintf (stderr, "BFS %s altered seed %d, %d->%d\n",
		 alg_type, i, seeds[i], inseeds[i]);

  identify_comm (num_communities, comm_root, comm_size, membership, member_flag);
}

static void
run_PR (int *num_communities, attr_id_t *membership, attr_id_t *member_flag)
{
  attr_id_t num_comm;
  double mod, cond, cc;
  attr_id_t inseeds[MAX_NUM_SEEDS];

  for (int i = 0; i < num_seeds; ++i) inseeds[i] = seeds[i];

  OMP("omp parallel for")
    for (int i = 0; i < g->n; ++i)
      member_flag[i] = -9999;

#define PR_ALPHA 0.0625
#define PR_EPS FLT_EPSILON

  pagerank_community (g, inseeds, num_seeds, member_flag, PR_ALPHA, PR_EPS);

  if (getenv ("VERBOSE"))
    for (int i = 0; i < num_seeds; ++i)
      if (inseeds[i] != seeds[i])
	fprintf (stderr, "PR %s altered seed %d, %d->%d\n",
		 alg_type, i, seeds[i], inseeds[i]);

  identify_comm (num_communities, comm_root, comm_size, membership, member_flag);
}

static void
run_RW (int *num_communities, attr_id_t *membership, attr_id_t *member_flag)
{
  attr_id_t inseeds[MAX_NUM_SEEDS];

  for (int i = 0; i < num_seeds; ++i) inseeds[i] = seeds[i];

  OMP("omp parallel for")
    for (int i = 0; i < g->n; ++i)
      member_flag[i] = -9999;

#define NSWEEP 60
  andersen_lang (g, inseeds, num_seeds, member_flag, NSWEEP);

  if (getenv ("VERBOSE"))
    for (int i = 0; i < num_seeds; ++i)
      if (inseeds[i] != seeds[i])
	fprintf (stderr, "RW %s altered seed %d, %d->%d\n",
		 alg_type, i, seeds[i], inseeds[i]);

  identify_comm (num_communities, comm_root, comm_size, membership, member_flag);
}

static int compare (const void * a, const void * b)
{
	return ( *(const int*)a - *(const int*)b );
}

int main (int argc, char* argv[])
{
	char *graph_type, *seed_file;
	int i;

	/* Initialize the random number generator */
	srand48 (time (NULL));

	/* Parse command line options  */
	if (argc < 3) {
	  if (argc > 1 && 0 == strcasecmp ("--header", argv[1])) {
	    dump_output (1);
	    return 0;
	  } else {
	    fprintf (stderr, "./eval_seed_community_detection "
		     "<seed_file> <alg_type> [<round_robin: %d> "
		     "<stop-on-merge: %d> <tolerance: %e> <threshold %d>].\n",
		     round_robin, stop_on_merge, tol, threshold);
	    return -1;
	  }
	}

	seed_file = argv[1];
	alg_type   = argv[2];

	if (argc > 3)
	  round_robin = atoi(argv[3]);
	if (argc > 4)
	  stop_on_merge = atoi(argv[4]);
	if (argc > 5)
	  tol = atof(argv[5]);
	if (argc > 6)
	  threshold = atoi(argv[6]);
	if(threshold == 0) /* Means no or inf threshold */
		threshold = INFTY;

	FILE * infile = fopen(seed_file, "r");
	infilename = (char *) malloc(sizeof(char)*255);

	fscanf(infile, "%s", infilename);

	fscanf(infile, "%d", &num_seeds);
	if (num_seeds > MAX_NUM_SEEDS) {
	  fprintf (stderr, "Requested number of seeds %d > max %d\n", num_seeds, MAX_NUM_SEEDS);
	  return -1;
	}
	for (int k = 0; k < num_seeds; ++k)
	  fscanf (infile, "%d", &seeds[k]);

	/* Generate graph */
	graph_type = (char *) calloc(500, sizeof(char));
	graph_ext_check(infilename, graph_type);
	g = (graph_t *) malloc(sizeof(graph_t));
	assert(g != NULL);
	graph_gen(g, infilename, graph_type);

OMP("omp parallel for")
	for(i=0; i<g->n; i++)
	{
		int start = g->numEdges[i];
		int end = g->numEdges[i+1];

		qsort(g->endV + start, end-start, sizeof(int), compare);
	}

	/* Run algorithm */
	membership = (attr_id_t *) malloc(g->n * sizeof(attr_id_t));
	member_flag = (attr_id_t *) malloc(g->n * sizeof(attr_id_t));

	if (0 == strcasecmp ("CNM", alg_type) || 0 == strcasecmp ("MB", alg_type))
	  run_MB_CNM (&num_communities, membership, member_flag);
	else if (0 == strcasecmp ("BFS", alg_type))
	  run_BFS (&num_communities, membership, member_flag);
	else if (0 == strcasecmp ("PR", alg_type))
	  run_PR (&num_communities, membership, member_flag);
	else if (0 == strcasecmp ("RW", alg_type))
	  run_RW (&num_communities, membership, member_flag);
	else {
	  fprintf (stderr, "Unrecognized algorithm %s\n", alg_type);
	  return -1;
	}

	clst = get_single_community_clustering_coefficient(g, member_flag, 1);
	modularity = get_single_community_modularity (g, member_flag, 1);
	conductance = get_single_community_conductance(g, member_flag, 1);

	dump_output (0);

	if (getenv ("DUMP_COMM")) {
	  FILE* cout = NULL;
	  cout = fopen (getenv ("DUMP_COMM"), "w+");
	  if (cout) {
	    for (i = 0; i < g->n; ++i)
	      fprintf (cout, "%ld\n", (long)membership[i]);
	    fclose (cout);
	  }
	}
	return 0;
}
