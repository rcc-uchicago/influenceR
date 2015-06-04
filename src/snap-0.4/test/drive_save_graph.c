#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <alloca.h>

#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "graph_partitioning.h"
#include "utils.h"

static int
comp_fn(const void *v1, const void *v2)
{
	const attr_id_t *p1, *p2;
	p1 = (attr_id_t *)v1;
	p2 = (attr_id_t *)v2;
	
	return (int)(*p1 - *p2);
}

int main (int argc, char* argv[])
{
  char *infilename, *outfilename;
  char *graph_type;
  graph_t *g;

  /* Parse command line options */
  if(argc < 2){
     fprintf (stderr, "Usage: %s <graph>\n", argv[0]);
     return -1;
  }

  infilename = argv[1];
  outfilename = "save.gr";

  /* Generate graph */
  graph_type = (char *) alloca (500 * sizeof(char));
  memset (graph_type, 0, 500);
  graph_ext_check(infilename, graph_type);
  g = (graph_t *) alloca (sizeof (graph_t));
  assert(g != NULL);
  graph_gen(g, infilename, graph_type);

  if (!g->undirected) {
    fprintf (stderr, "Only undirected graphs implemented.");
    return -1;
  }

  /*dump the graph to file */
  FILE *fp = fopen (outfilename, "w");
  save_undir_unwgt_graph(fp, g);
  fclose(fp);
  fprintf (stderr, "Saved the graph to file '%s'...\n", outfilename);

  /* Test whether we have the same graph as we wrote or not */
  attr_id_t i, j, n, start_iter, end_iter, degree;
  attr_id_t *tmp, *tmp2;
  graph_t *sg;
  char *sgraph_type;

  sgraph_type = (char*) alloca (500 * sizeof(char));
  memset (sgraph_type, 0, 500);
  graph_ext_check(outfilename, sgraph_type);
  sg = (graph_t *) alloca (sizeof(graph_t));
  assert(sg != NULL);
  graph_gen(sg, outfilename, sgraph_type);
  fprintf(stderr, "Read in the saved graph...\n");

  /* Compare the number of vertices and the number of edges */
  assert (g->n == sg->n);
  assert (g->m == sg->m);

  tmp = (attr_id_t*) malloc (g->m * sizeof(attr_id_t));
  tmp2 = (attr_id_t*) malloc (g->m * sizeof(attr_id_t));
  for( i=0; i<n; i++){
	  assert ( (g->numEdges[i+1] - g->numEdges[i]) == (sg->numEdges[i+1] - sg->numEdges[i]));
      start_iter = g->numEdges[i];
	  end_iter = g->numEdges[i+1];
	  degree = end_iter - start_iter;

	 memcpy(tmp, g->endV+start_iter, degree * sizeof(attr_id_t));
	 qsort(tmp, degree, sizeof(attr_id_t), comp_fn);
	 memcpy(tmp2, sg->endV+start_iter, degree * sizeof(attr_id_t));
	 qsort(tmp2, degree, sizeof(attr_id_t), comp_fn);
     for(j=0; j<degree; j++){
		assert(tmp[j] == tmp2[j]); 
	 }
  }

  free (tmp);
  free (tmp2);
  free_graph(g);
  free_graph(sg);

  return 0;
}
