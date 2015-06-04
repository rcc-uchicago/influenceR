#include "graph_defs.h"
#include "graph_gen.h"
#include "graph_kernels.h"
#include "graph_metrics.h"
#include "graph_partitioning.h"
#include "utils.h"

#define DEBUG

int main(int argc, char** argv) {
    char *infilename, *outfilename, *graph_type, *alg_type;
    FILE* fp;
    graph_t* g;
    int curArgIndex;
    attr_id_t *membership;
    int num_communities;
    double modularity, mod_val, time0;
    long i;
    int run_approxBC;
    int digits;
    char *forstr;

    /* Step 1: Parse command line arguments */
    if (argc < 3) {
        fprintf(stdout, "\nUsage: %s -infile <graph filename>"
                " -graph <graph type> -outfile <output filename>)"
                " -alg <algorithm type>\n\n",
                "eval_modularity_greedy_agglomerative");
        fprintf(stdout, "Algorithm type can be one of the following:\n"
                "CNM         -- greedy agglomerative strategy of "
                "Clauset, Newman and Moore.\n"
                "MB          -- McCloskey and Bader: "
                "normalize CNM Qij's by stddev(Qij).\n"
                "RAT         -- normalize CNM by size ratios.\n"
                "MBRAT       -- normalize McCloskey-Bader by size ratios.\n"
                "LIN         -- McCloskey and Bader: weighted graphs "
                "with linear model.\n"
#if 0
                "WT1/WT2/WT3 -- Wakita and Tsurami's"
                "consolidation ratio heuristics.\n"
                "DDA         -- Normalized modularity heuristics "
                "by Danon, Diaz-Guilera, and Arenas.\n"
                "CCA         -- An agglomerative clustering heuristic "
                "based on local clustering coefficient.\n"
#endif
                "\n");

        usage_graph_options();
        exit(-1);
    }

    curArgIndex = 0;
    run_approxBC = 0;
    infilename = (char *) calloc(500, sizeof(char));
    outfilename = (char *) calloc(500, sizeof(char));
    graph_type = (char *) calloc(500, sizeof(char));
    alg_type = (char *) calloc(500, sizeof(char));

    strcpy(alg_type,"CNM");
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

        if (strcmp(argv[curArgIndex],"-alg")==0) {
#ifdef DEBUGXX
            printf("argv: %s (len: %d)\n",argv[curArgIndex+1], strlen(argv[curArgIndex+1]));
#endif
            strcpy(alg_type, argv[++curArgIndex]);
            printf("alg_type: %s\n",alg_type);
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

    membership = (attr_id_t *) malloc(g->n * sizeof(attr_id_t));
    assert(membership != NULL);
    modularity = 0.0;
    num_communities = 0;
    time0 = get_seconds();
    modularity_greedy_agglomerative(g, alg_type, membership, 
            &num_communities, &modularity);
    time0 = get_seconds() - time0;

    mod_val = get_community_modularity(g, membership, num_communities);

    /* Step 4: Write output to file */
    fp = fopen(outfilename, "w");
    fprintf(fp, "Number of communities: %d\n", num_communities);
    fprintf(fp, "Modularity score: %f (full), %f (w/o dup)\n", mod_val,
            modularity);
    fprintf(fp, "\n<Vertex ID> <Community ID>\n\n");

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
                fprintf(fp, forstr, i, membership[i]);  
            else
                fprintf(fp, forstr, i+1, membership[i]);  
        }

    }
    fprintf(stderr, "Modularity: %f (full), %f (w/o dup)\n", mod_val,
            modularity);
    fprintf(stderr, "Number of communities: %d\n", num_communities);
    fprintf(stderr, "Execution Time: %f seconds\n", time0);

#if 0
    {
        int *marked;
        long i, j, v;
        long edgecount, eIdx, idx, currLabelCount, currLabelSources;

        typedef struct {
            int from;
            int from_label;
            int to;
            int to_label;
        } cross_edge_t;

        cross_edge_t *edges;

        int cross_edge_compar(const void *v1, const void *v2) {
            const cross_edge_t *p1, *p2;
            p1 = (cross_edge_t *)v1;
            p2 = (cross_edge_t *)v2;

            if ((p1->from == p2->from) && (p1->from_label == p2->from_label) && 
                    (p1->to   == p2->to)   && (p1->to_label   == p2->to_label))
                return 0;

            if (p1->from_label < p2->from_label)
                return -1;

            if (p1->from_label > p2->from_label)
                return +1;

            return (p1->from < p2->from ? -1 : +1);
        }


        marked = (int *)malloc(g->n * sizeof(int));
        assert(marked != NULL);

        for (i=0; i<g->n; i++)
            marked[i] = (g->numEdges[i+1] > g->numEdges[i] ? 1:  0); /* Degree > 0 */

        edgecount = 0;
        for (i=0; i<g->n; i++) {
            if (marked[i]) {
                for (j=g->numEdges[i] ; j<g->numEdges[i+1] ; j++) {
                    v = g->endV[j];
                    if (!marked[v]) fprintf(stderr,"ERROR: v %d not marked\n",v);
                    if (membership[i] != membership[v]) {
#if 0
                        fprintf(stderr,"TEST %4d (%4d) != %4d (%4d)\n",
                                i, membership[i], v, membership[v]);
#endif
                        edgecount++;
                    }
                }
            }
        }

        if (edgecount > 0) {

            edges = (cross_edge_t *)malloc(edgecount * sizeof(cross_edge_t));
            assert(edges != NULL);


            eIdx = 0;
            for (i=0; i<g->n; i++) {
                if (marked[i]) {
                    for (j=g->numEdges[i] ; j<g->numEdges[i+1] ; j++) {
                        v = g->endV[j];
                        if (membership[i] != membership[v]) {
                            edges[eIdx].from = i;
                            edges[eIdx].from_label = membership[i];
                            edges[eIdx].to   = v;
                            edges[eIdx].to_label   = membership[v];
                            eIdx++;
                        }
                    }
                }
            }

            qsort(edges, edgecount, sizeof(cross_edge_t), cross_edge_compar);

#if 1
            for (i=0 ; i<edgecount ; i++)
                fprintf(stderr,"crossedge[%4d]: %4d %4d %4d %4d\n",
                        i, edges[i].from, edges[i].from_label, edges[i].to, edges[i].to_label);
#endif

            fprintf(stdout,"\n");

            idx = 0;
            currLabelCount   = 1;
            currLabelSources = 1;
            for (i=1 ; i<edgecount ; i++) {
                if (edges[i].from_label == edges[idx].from_label) {
                    currLabelCount++;
                    if (edges[i].from != edges[i-1].from)
                        currLabelSources++;
                }
                else {
                    if (currLabelCount > 2) {
#ifdef DEBUG
                        fprintf(stderr,"summary: label: %4d  count: %4d  sources: %4d  FILE: %4d\n",
                                edges[i-1].from_label, currLabelCount, currLabelSources, edges[i-1].from);
#else
                        if (currLabelSources==1)
                            fprintf(stdout,"T NODE: %4d\n", edges[i-1].from);
#endif
                    }
                    idx = i;
                    currLabelCount = 1;
                    currLabelSources = 1;
                }
            }
            if ((idx <edgecount) && (currLabelCount > 2)) {
#ifdef DEBUG
                fprintf(stderr,"summary: label: %4d  count: %4d  sources: %4d  FILE: %4d\n",
                        edges[idx].from_label, currLabelCount, currLabelSources, edges[idx].from);
#else
                if (currLabelSources==1)
                    fprintf(stdout,"T NODE: %4d\n", edges[idx].from);
#endif
            }      


            free(edges);

        }
        else {
            fprintf(stdout,"No edges\n");
        }

        free(marked);
    }
#endif

    fclose(fp);

    free(forstr);

    /* Step 5: Clean up */
    free(infilename);
    free(outfilename);
    free(graph_type);
    free(alg_type);
    free_graph(g);
    free(g);
    free(membership);

    return 0;
}
