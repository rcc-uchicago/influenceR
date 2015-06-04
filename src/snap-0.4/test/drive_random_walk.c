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

int main (int argc, char* argv[])
{
     char *infilename;
     char *graph_type;
     graph_t *g;
     attr_id_t *seeds;
     unsigned long i, nseeds = 3, walklen = 3;

     /* Initialize the random number generator */
     srand48 (time (NULL));

     /* Parse command line options  */
     if (argc > 4) {
	  fprintf (stderr, "%s "
		   "<graph> <nseeds = %lu> <walk-len = %lu>\n",
		   argv[0], nseeds, walklen);
	  return -1;
     }
    
     infilename = argv[1];
     if (argc > 2) {
	  nseeds = strtoul (argv[2], NULL, 10);
	  if (!nseeds /*|| nseeds == ULLONG_MAX*/) {
	       perror ("Failed to parse nseeds.");
	       return -1;
	  }
     }
     if (argc > 3) {
	  walklen = strtoul (argv[3], NULL, 10);
	  if (!walklen /*|| walklen == ULLONG_MAX*/) {
	       perror ("Failed to parse walk-len.");
	       return -1;
	  }
     }

     /* Generate graph */
     graph_type = (char *) alloca (500 * sizeof (char));
     memset (graph_type, 0, 500);
     graph_ext_check(infilename, graph_type);
     g = (graph_t *) alloca (sizeof (graph_t));
     assert(g != NULL);
     graph_gen(g, infilename, graph_type);

     seeds = alloca (nseeds * sizeof (*seeds));
     generate_random_walk_seeds(g, nseeds, walklen, seeds);

     printf ("%s\n%lu\n%lu", infilename, nseeds, (unsigned long)seeds[0]);
     for (i = 1; i < nseeds; ++i)
	  printf (" %lu", (unsigned long)seeds[i]);
     printf ("\n");

     return 0;
}
