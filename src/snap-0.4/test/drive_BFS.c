#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "utils.h"

int main(int argc, char** argv) {

    char *infilename, *outfilename, *graph_type;
    FILE* fp;
    graph_t* g;

    long src;
    int curArgIndex;
    int est_diameter;
    long num_vertices_visited;
    
    int proc_pid;
    
    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s (-src <vertex ID (0 to n-1)>) -infile <graph filename>"
                " (-graph <graph type> -outfile <output filename>)\n\n", "eval_BFS");
        
        usage_graph_options();
        
        fprintf(stdout, "Using the -src option, specify the vertex ID to run BFS from. A random vertex is selected if the src is not specified.\n\n");
        exit(-1);
    }

    curArgIndex = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));

    proc_pid = getpid();
    sprintf(outfilename, "results.%d.txt", proc_pid);

    src = -1;

    while (curArgIndex < argc) {
        
        if (strcmp(argv[curArgIndex],"-src")==0) {
            src = atol(argv[++curArgIndex]);
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

    print_snap_header(stdout);
    
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

    print_snap_header(fp);

    graph_ext_check(infilename, graph_type);

    fprintf(stdout, "Input graph file   : %s\n", infilename);
    fprintf(stdout, "Output file  : %s\n", outfilename);

    fprintf(fp, "Input graph file    : %s\n", infilename);
    fprintf(fp, "Output file   : %s\n", outfilename);
    
    /* Step 2: Generate graph */
    g = (graph_t *) malloc(sizeof(graph_t));
    graph_gen(g, infilename, graph_type);
   
    print_graph_header(stdout, g, "BFS");
    print_graph_header(fp, g, "BFS");
    
    if (src == -1)
       src = lrand48() % g->n;

    assert((src >= 0) && (src < g->n));

    fprintf(stdout, "  Source vertex      : %ld\n\n", src);
    fprintf(fp    , "  Source vertex      : %ld\n\n", src);
    
    /* Step 3: Run algorithm */

    /* Assuming a low diameter graph */
    /* change the est_diameter value for high diameter graphs */
    est_diameter = 100;
    num_vertices_visited = BFS_parallel_frontier_expansion(g, src, est_diameter);

    /* Step 4: Write output to file */
    fprintf(fp, "  Breadth-first search from vertex %ld\n", src); 
    fprintf(fp, "  Number of vertices visited: %ld\n\n", num_vertices_visited);
    fprintf(stdout, "  Breadth-first search from vertex %ld\n", src); 
    fprintf(stdout, "  Number of vertices visited: %ld\n\n", 
            num_vertices_visited);
    
    /* Step 5: Clean up */
    free(infilename);
    free(outfilename);
    free(graph_type);

    free_graph(g);
    free(g);
    fclose(fp);
    return 0;
}
