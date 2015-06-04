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
    
    attr_id_t* component_num;
    attr_id_t i, bcc_num;
    int digits;
    char *forstr;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                " (-graph <graph type> -outfile <output filename>)\n\n", "eval_biconnected_components");
        
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
        
    /* Step 3: Run algorithm */ 

    component_num = (attr_id_t *) malloc(g->m*sizeof(attr_id_t));
    bcc_num = biconnected_components(g, component_num);

    /* Write output to file */
    fp = fopen(outfilename, "w");
    fprintf(fp, "Number of biconnected components: %d\n", bcc_num);
    fprintf(fp, "\n<Vertex ID> <Biconnected Component ID>\n\n");

    if (g->n > 1)
      digits =  1 + (int)floor(log ((double)g->n) / log (10.0));
    else
      digits = 1;

    forstr = (char *)malloc(256*sizeof(char));
    assert(forstr != NULL);

    assert(digits < 32);

    sprintf(forstr, "%%%dld %%%dd\n", digits, digits);

    for (i=0; i<g->n; i++) {
      if (g->numEdges[i+1] - g->numEdges[i] > 0) {
        if (g->zero_indexed)
	  fprintf(fp, forstr, i, component_num[i]);  
        else
	  fprintf(fp, forstr, i+1, component_num[i]);  
      }

    }
    fflush(fp);
    fclose(fp);
    free(forstr);

 
    /* Step 4: Clean up */
    free(component_num);

    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);

    return 0;
}
