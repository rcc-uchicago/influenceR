#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "graph_partitioning.h"
#include "utils.h"

int main(int argc, char** argv) {

    char *infilename, *outfilename, *graph_type;
    FILE* fp;
    graph_t* g;
    int curArgIndex;
    attr_id_t *membership;
    int num_communities;
    double modularity;
    int with_klin;
    long i;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                " (-klin -graph <graph type> -outfile <output filename>)\n\n",
                "eval_modularity_spectral");
        
        usage_graph_options();
        exit(-1);
    }

    curArgIndex = 0;
    with_klin = 0;
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
 
        if (strcmp(argv[curArgIndex],"-klin")==0) {
            with_klin = 1;
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
    
    /* Fix in future version of code */    
    if (g->undirected) {
	    g->m /=2;
    }
    
    membership = (attr_id_t *) malloc(g->n * sizeof(attr_id_t));
    assert(membership != NULL);

    if (with_klin) { 
        fprintf(stderr, "Running the spectral algorithm for modularity "
            "optimization (with Kernighan-Lin refinement)\n");
        modularity_spectral(g, membership, &num_communities, 1);
	    computeModularityValue(g, membership, num_communities, &modularity);
    } else {
        fprintf(stderr, "Running the spectral algorithm for modularity "
            "optimization (without Kernighan-Lin refinement)\n");
        modularity_spectral_wo_klin(g, membership,&num_communities);
	    computeModularityValue(g, membership, num_communities, &modularity);
    }

    /* Step 4: Write output to file */
    fp = fopen(outfilename, "w");
    fprintf(fp, "Number of communities: %d\n", num_communities);
    fprintf(fp, "Modularity score: %f\n", modularity);
    fprintf(fp, "\n<Vertex ID> <Community ID>\n\n");
    for (i=0; i<g->n; i++) {
        if (g->zero_indexed)
            fprintf(fp, "%ld %d\n", i, membership[i]);  
        else
            fprintf(fp, "%ld %d\n", i+1, membership[i]);  

    }
    fclose(fp);


    /* Step 5: Clean up */
    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);
    free(membership);

    return 0;
}
