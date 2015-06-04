#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"

int main(int argc, char** argv) {

    char *infilename, *outfilename, *graph_type;
    FILE* fp;
    graph_t* g;

    int curArgIndex;
    long i;
    int vc_weighted_size, vc_unweighted_size;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                " (-graph <graph type> -outfile <output filename>)\n\n",
                "eval_vertex_cover");
        
        usage_graph_options();
        exit(-1);
    }

    curArgIndex = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));

    strcpy(outfilename, "output.txt");

    while (curArgIndex < argc) {
        
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

    /* Step 2: Generate graph */
    g = (graph_t *) malloc(sizeof(graph_t));
    graph_gen(g, infilename, graph_type);
   
    fprintf(stdout, "Number of vertices     : %ld\n", g->n);
    if (g->undirected)
        fprintf(stdout, "Number of edges        : %ld\n\n", g->m/2);
    else 
        fprintf(stdout, "Number of edges        : %ld\n\n", g->m);

    if (g->undirected == 0) {
        fprintf(stderr, "Error: the graph has to be undirected.\n");
        fprintf(stderr, "Please check input file.\n");
        exit(-1);
    }

    /* Step 3: Run algorithm */
	g->dbl_weight_v = (double *) malloc(g->n * sizeof(double));
	
    if (g->undirected) {
	    g->m /=2;
    }

	for (i=0; i<g->n; i++) {
		g->dbl_weight_v[i] = 1;
	}
	
    vc_weighted_size = vertex_cover_weighted(g);
	
	vc_unweighted_size = vertex_cover_unweighted(g);

    /* Step 4: Write output to file */
    fp = fopen(outfilename, "w");
    fprintf(fp, "Vertex cover size (weighted): %d\n", vc_weighted_size);
    fprintf(fp, "Vertex cover size (unweighted): %d\n", vc_unweighted_size);
    
    /* Step 5: Clean up */
    free(g->dbl_weight_v);
    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);

    return 0;
}
