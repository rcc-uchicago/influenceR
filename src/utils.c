#include <time.h>
#include <unistd.h>

#include "utils.h"

double snap_runtime;

double get_seconds() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return (double) (tp.tv_sec + ((1e-6)*tp.tv_usec));
}

void usage_graph_options() {
    fprintf(stdout, "The input file name is a required argument. We try to determine the graph");
    fprintf(stdout, " format from the file extension; optionally, you can specify the graph type with");
    fprintf(stdout, " the -graph option.\n\n");
    fprintf(stdout, "Supported graph file formats (for use with the -graph option):\n");
    fprintf(stdout, "snap    (.gr)       SNAP file format (default file format supported by this package).\n");
    fprintf(stdout, "dimacs  (.dim)      graph format used in the 9th DIMACS Shortest Paths Challenge.\n");
    fprintf(stdout, "GML     (.gml)      a limited implementation of the GML graph format (only supports undirected graphs, edge weights represented as doubles.)\n");
#if 0
    /* Commenting out unimplemented options */
    fprintf(stdout, "metis   (.met)      graph representation used by the Metis partitioning package.\n");
    fprintf(stdout, "graphml (.graphml)  GraphML representation.\n");
    fprintf(stdout, "sqmesh  (.sqm)      a synthetic 2D square mesh generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "lmesh   (.lm)       a synthetic 2D long mesh generator is invoked."
            " The input file specifies the generator parameters.\n");
#endif    
    fprintf(stdout, "rand    (.rnd)      a synthetic random graph generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "rmat    (.rmat)     a synthetic scale-free graph generator is invoked."
            " The input file specifies the generator parameters.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "The default output filename is \"results.PID.txt\"."
            " (PID is the process ID). "
            "Specify an alternate file name with the -outfile option.\n");
}

void graph_ext_check(char* filename, char* graph_type) {

    char *ptr, *tmp; 

    ptr = filename;

    /* Find the last "." in the file name */
    do {
        tmp = strstr(ptr+1, ".");
        if (tmp != NULL) 
            ptr = tmp;
    }  while (tmp != NULL);
    ptr++;   

    if (graph_type[0] == 0) {
        if (strcmp(ptr, "gr") == 0) {
            strcpy(graph_type, "snap");
        }

        if (strcmp(ptr, "dim") == 0) {
            strcpy(graph_type, "dimacs");    
        }

        if (strcmp(ptr, "gml") == 0) {
            strcpy(graph_type, "gml");    
        }

        if (strcmp(ptr, "rnd") == 0) {
            strcpy(graph_type, "rand");    
        }

        if (strcmp(ptr, "rmat") == 0) {
            strcpy(graph_type, "rmat");    
        }
#if 0
        if (strcmp(ptr, "met") == 0) {
            strcpy(graph_type, "metis");    
        }

        if (strcmp(ptr, "graphml") == 0) {
            strcpy(graph_type, "graphml");    
        }

        if (strcmp(ptr, "sqmesh") == 0) {
            strcpy(graph_type, "sqmesh");    
        }

        if (strcmp(ptr, "lmesh") == 0) {
            strcpy(graph_type, "lmesh");    
        }
#endif
    }

}

void prefix_sums(attr_id_t *input, attr_id_t* result, attr_id_t* p, attr_id_t n) {

#ifdef _OPENMP
    attr_id_t i, j, r, start, end, add_value;
    int tid, nthreads;

    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    r = n/nthreads;

    result[0] = 0;

OMP("omp for")
    for (i=1; i<n+1; i++)
        result[i] = input[i-1];

    start =  tid*r + 1;
    end   = (tid+1)*r;

    if (tid == nthreads-1)
        end = n+1;

    for (j=start; j<end; j++)
        result[j] = input[j-1] + result[j-1];

    p[tid] = result[end-1];

OMP("omp barrier")

    if (tid == 0) {
        for (j=1; j<nthreads; j++)
            p[j] += p[j-1];
    }

OMP("omp barrier")

    if (tid>0) {
        add_value=p[tid-1];
        for (j=start-1; j<end; j++)
            result[j] += add_value;
    }

OMP("omp barrier")
#else

    attr_id_t i;
    result[0] = 0;
    for (i=1; i<n+1; i++) {
        result[i] = result[i-1] + input[i-1];
    }
#endif

}



void print_graph(graph_t* G)
{
    attr_id_t i,j;
    attr_id_t degree,start;

    fprintf(stdout, "\n Printing graph to stdout\n");

    fprintf(stdout, "Number of vertices: %ld, edges: %ld\n", G->n, G->m);
    if(G->undirected) 
        fprintf(stdout, "Graph is undirected\n");
    else
        fprintf(stdout, "Graph is directed\n");

    fprintf(stdout, "Graph weight type is %d\n",G->weight_type);

#if ENABLE_64BIT_VID

    for(i = 0; i < G->n; i++) {
        degree = G->numEdges[i+1] - G->numEdges[i];
        start = G->numEdges[i];
        fprintf(stdout, "Degree of vertex i %ld is %ld." 
                "Its neighbours are:", i, degree);

        for(j=0;j<degree;j++) {
            if(G->weight_type == 1)
                fprintf(stdout, "%ld [%d], ",G->endV[start+j], 
                        G->int_weight_e[start+j]);
            else if(G->weight_type == 2) 
                fprintf(stdout, "%ld [%ld], ",G->endV[start+j], 
                        G->l_weight_e[start+j]);
            else if(G->weight_type == 3) 
                fprintf(stdout, "%ld [%f], ",G->endV[start+j], 
                        G->fl_weight_e[start+j]);
            else if(G->weight_type == 4) 
                fprintf(stdout, "%ld [%f], ",G->endV[start+j], 
                        G->dbl_weight_e[start+j]);
        }
        fprintf(stdout, "\n");
    }

#else

    for(i = 0; i < G->n; i++) {
        degree = G->numEdges[i+1] - G->numEdges[i];
        start = G->numEdges[i];
        fprintf(stdout, "Degree of vertex i %d is %d." 
                "Its neighbours are:", i, degree);

        for(j=0;j<degree;j++) {
            if(G->weight_type == 1)
                fprintf(stdout, "%d [%d], ",G->endV[start+j], 
                        G->int_weight_e[start+j]);
            else if(G->weight_type == 2) 
                fprintf(stdout, "%d [%ld], ",G->endV[start+j], 
                        G->l_weight_e[start+j]);
            else if(G->weight_type == 3) 
                fprintf(stdout, "%d [%f], ",G->endV[start+j], 
                        G->fl_weight_e[start+j]);
            else if(G->weight_type == 4) 
                fprintf(stdout, "%d [%f], ",G->endV[start+j], 
                        G->dbl_weight_e[start+j]);
        }
        fprintf(stdout, "\n");
    }
#endif
}

void print_snap_header(FILE *outfile) {
    fprintf(outfile, "\n"
            "************************************************\n");
    fprintf(outfile, "  SNAP Complex Network Analysis Framework v0.3\n");
    fprintf(outfile, "  Authors: Kamesh Madduri, David A. Bader  \n");
    fprintf(outfile, "  Last Updated: March 2009\n");
    fprintf(outfile, "  http://snap-graph.sourceforge.net\n");
    fprintf(outfile, "************************************************\n");
}


void print_graph_header(FILE *outfile, graph_t *g, const char *problem_type) {
    fprintf(outfile, "  Number of vertices : %ld\n", g->n);
    if (g->undirected) {
        fprintf(outfile, "  Number of edges    : %ld\n", g->m/2);
        fprintf(outfile, "  Undirected network.\n");
    } else {
        fprintf(outfile, "  Number of edges    : %ld\n", g->m);
        fprintf(outfile, "  Directed network.\n");
    }
    /*
       if (g->weight_type == 0)
       fprintf(outfile, "  Unweighted network.\n");
       else if (g->weight_type == 1)
       fprintf(outfile, "  32-bit integer weight.\n");
       else if (g->weight_type == 2)
       fprintf(outfile, "  Long integer (%d bytes) weights.\n", sizeof(long));
       else if (g->weight_type == 3)
       fprintf(outfile, "  Single-precision float weights.\n");
       else if (g->weight_type == 4)
       fprintf(outfile, "  Double-precision float weights.\n");
     */ 
    fprintf(outfile, "------------------------------------------------\n");
    fprintf(outfile, "  Problem type       : %s\n", problem_type); 
}

static int
comp_fn(const void *v1, const void *v2)
{
     return (int)(*(attr_id_t*)v1 - *(attr_id_t*)v2);
}

/** Saves the graph in snap(.gr) format. Assumes graph to be
  undirected, unweighted and zero indexed
*/
void
save_undir_unwgt_graph(FILE *outfile, graph_t *g)
{
     attr_id_t i, j, v, prev_v, start_iter, end_iter, m, n, degree, count;
     attr_id_t *endV, *numEdges, *tmp, *tmp2;

     n = g->n;
     m = g->m/2;
     endV = g->endV;
     numEdges = g->numEdges;
     tmp = (attr_id_t*) malloc (g->m * sizeof(attr_id_t));

     /* Assuming that the graph contains vertices labelled starting from 0 */
     assert (outfile);
     fprintf (outfile, "p %ld %ld u u 0\n", n, m);
     for(i=0; i<n; i++){
	  start_iter = numEdges[i];
	  end_iter = numEdges[i+1];
	  degree = end_iter - start_iter;
	  prev_v = -1; count = 0;
	  memcpy(tmp, endV+start_iter, degree * sizeof(attr_id_t));
	  qsort(tmp, degree, sizeof(attr_id_t), comp_fn);
	  for(j=0; j<degree; j++){
	       v = tmp[j];
	       if(prev_v == v)
		    count++;
	       else
		    count = 1;
	       prev_v = v;

	       if(i < v){
		    fprintf (outfile, "%ld %ld\n", i, v);
	       }
	       else if (i == v){
		    if(count % 2 == 1)
			 fprintf (outfile, "%ld %ld\n", i, v);
	       }
	  }
     }
}

static attr_id_t
random_walk(graph_t *g, attr_id_t V, int num_levels)
{
     int i;
     attr_id_t index, deg_V, rand_num;

     for(i=0; i<num_levels; i++) {
	  // number of vertices adjacent to the vertex V i.e. degree[V]
	  deg_V = g->numEdges[V + 1] - g->numEdges[V];
	  if(deg_V == 0) {
	       //if degree of V is zero then new random vertex will be taken as the source	       return -2;
	  }
	  rand_num = (attr_id_t) (drand48() * deg_V);
	  index = g->numEdges[V] + rand_num;
	  V = g->endV[index];
     }

     return(V);

}

/** Generate a seed set using random walks.

    @param g Graph
    @param num_seeds Number of seeds to generate
    @param num_levels Number of steps in the walk
    @param seeds Output seeds
*/
void
generate_random_walk_seeds(graph_t *g, int num_seeds, int num_levels, attr_id_t *seeds)
{
     attr_id_t S;
     int size_seeds = 0;
     int flag = 0;
     int i,j;
     int num_repeats = 0;
     attr_id_t random_vertex;
     time_t t0;

     /* Initialize the random number generator */
     t0 = ((int)time (NULL)) ^ ((int)getpid ());
     srand48 (t0);

     //Generating a random vertex to begin the random walk
     random_vertex = (attr_id_t) (drand48() * g->n);

     for(i=0; i<num_seeds; i++)
	  seeds[i] = -1;

     /*printf("\nInitialising the Random Walk within the graph========\n");*/
     while(size_seeds < num_seeds) {
	  S = random_walk(g, random_vertex, num_levels);

	  if(S == -2) {
	       random_vertex = (attr_id_t) (drand48() * g->n);
	       continue;
	  }

	  for(j=0; j<size_seeds; j++)
	       if(S == seeds[j])
		    flag = 1;

	  if(flag == 0) {
	       size_seeds++;
	       seeds[size_seeds-1] = S;
	  } else {
	       flag = 0;
	       num_repeats++;
	       if(num_repeats > 100) {
		    //if the same vertex gets repeated many times then choose another source vertex
		    /*printf("\n Number of repeats are too large...Generating another random vertex\n");*/
		    random_vertex = (attr_id_t) (drand48() * g->n);
		    /*printf("\n The new random vertex is %d \n", (int)random_vertex);*/
	       }
	  }
     }
}
