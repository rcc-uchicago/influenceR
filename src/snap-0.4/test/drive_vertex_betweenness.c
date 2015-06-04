#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"

typedef struct {
  double val;
  long   idx;
} vBC_t;

int compare_vBC(const void *v1, const void *v2) {

  const vBC_t *p1, *p2;
  p1 = (vBC_t *)v1;
  p2 = (vBC_t *)v2;

  if (p1->val == p2->val) {
    if (p1->idx == p2->idx) 
      return 0;
    else 
      if (p1->idx < p2->idx) 
	return -1;
      else
	return 1;
  }

  if (p1->val > p2->val) 
    return -1;
  else
    return 1;
}


int main(int argc, char** argv) {

    char *infilename, *outfilename, *graph_type;
    FILE* fp;    
    graph_t* g;
    double* vBC;

    vBC_t *vBetArray;
    

    int run_approx_BC;
    double sampling_val;
    int curArgIndex;
    long numSrcs;

    long i;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s (-approx <sampling percentage>) -infile <graph filename>"
                " (-graph <graph type> -outfile <output filename>)\n\n", "eval_vertex_betweenness");
        
        usage_graph_options();
        
        fprintf(stdout, "\nTo compute approximate betweenness centrality, use the -approx flag with the"
               " percentage of vertices to run graph traversals from. 100, for instance, corresponds to"
               " exact betweenness computation and 0.1 runs BFSes from a random sample of n/1000 vertices.\n\n");
        exit(-1);
    }

    curArgIndex = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));

    sampling_val = 100.0;
    strcpy(outfilename, "output.txt");

    run_approx_BC = 0;

    while (curArgIndex < argc) {
        
        if (strcmp(argv[curArgIndex],"-approx")==0) {
            run_approx_BC = 1;
            sampling_val = atof(argv[++curArgIndex]);
        }

        if (strcmp(argv[curArgIndex],"-infile")==0) {
            strcpy(infilename, argv[++curArgIndex]);
        }

        if (strcmp(argv[curArgIndex], "-outfile")==0) {
            strcpy(outfilename, argv[++curArgIndex]);
        } 
 
        if (strcmp(argv[curArgIndex], "-graph")==0) {
            strcpy(graph_type, argv[++curArgIndex]);
        } 
        curArgIndex++; 
    }

    fp = fopen(infilename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not open input file. Exiting ...\n");
        exit(-1);
    }
    fclose(fp);

    fp = fopen(outfilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error! Could not write to output file. Exiting ...\n");   
        exit(-1);
    }
    fclose(fp);

    graph_ext_check(infilename, graph_type);

    fprintf(stdout, "\n");
    fprintf(stdout, "Input Graph File    : %s\n", infilename);
    fprintf(stdout, "Output Graph File   : %s\n\n", outfilename);

    if ((sampling_val <= 0) || (sampling_val > 100)) {
       fprintf(stderr, "Error! Sampling percentage should be between 0 and 100. Exiting ...\n"); 
       exit(-1);
    }
    fprintf(stdout, "Sampling percentage : %f\n\n", sampling_val);


    /* Step 2: Generate graph */
    g = (graph_t *) malloc(sizeof(graph_t));
    graph_gen(g, infilename, graph_type);
   
    fprintf(stdout, "Number of vertices     : %ld\n", g->n);
    if (g->undirected)
        fprintf(stdout, "Number of edges        : %ld\n\n", g->m/2);
    else
        fprintf(stdout, "Number of edges        : %ld\n\n", g->m);
        
    /* Step 3: Run algorithm */
    vBC = (double *) calloc(g->n, sizeof(double));
    numSrcs = g->n * (sampling_val/100.0);
    if (numSrcs == 0) 
        numSrcs = 1;
    vertex_betweenness_centrality(g, vBC, numSrcs);
  
    /* Step 4: Write results to output file */ 
    fp = fopen(outfilename, "w"); 
    fprintf(fp, "Input file: %s\n", infilename);
    if (g->undirected)
        fprintf(fp, "n: %ld, m: %ld\n", g->n, g->m/2);
    else 
        fprintf(fp, "n: %ld, m: %ld\n", g->n, g->m);
    fprintf(fp, "numSrcs: %ld\n", numSrcs);
    fprintf(fp, "\n<BC score> <Vertex ID>\n\n");

    vBetArray = (vBC_t *)malloc(g->n * sizeof(vBC_t));
    assert(vBetArray != NULL);

    for (i=0; i<g->n; i++) {
      vBetArray[i].val = vBC[i];
      if (g->zero_indexed) {
	fprintf(fp, "%15.7f %ld\n", vBC[i], i);
	vBetArray[i].idx = i;
	
      }
      else {
	fprintf(fp, "%15.7f %ld\n", vBC[i], i+1);  
	vBetArray[i].idx = i+1;
      }
    }

    qsort(vBetArray, g->n, sizeof(vBC_t), compare_vBC);

    for (i=0; i<g->n; i++) {
      fprintf(fp, "BC %ld %15.7f\n", vBetArray[i].idx, vBetArray[i].val);
    }

    fclose(fp);


    free (vBetArray);

    /* Step 5: Clean up */
    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);
    free(vBC);

    return 0;
}
