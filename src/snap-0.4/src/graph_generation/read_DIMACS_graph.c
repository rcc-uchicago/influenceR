#include "graph_defs.h"
#include "graph_gen.h"

void read_DIMACS_graph(graph_t* G, char* filename) {

    char *buf;
    FILE *fp;

    long i;
    long n, m;
    long count, offset;
    int int_wt, *int_weight;
    long u, v;
    attr_id_t *src;
    attr_id_t *dest;
    attr_id_t *degree;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: Cannot open input file. Exiting ...\n");
        exit(1);
    }

    /* assuming no line is longer than 500 characters */ 
    buf = (char *) malloc(500*sizeof(char));
    if (buf == NULL) {
        fprintf(stderr, "Error: Could not allocate memory"
                " in file read routine."
                " Exiting ...\n");
    }

    /* DIMACS graphs are undirected, weighted, and 1-indexed */
    /* Edge weights are signed integers */

    while (fgets(buf, 500, fp) != NULL)  {

        /* c is the comment symbol */
        if (buf[0] == 'c')
            continue;

        if (buf[0] == 'p') {
            sscanf(buf, "%*c %*2s %ld %ld", &n, &m);
        }
#if 1
        fprintf(stderr,"n: %d m: %d\n",n,m);
        fflush(stderr);
#endif
        assert(n>0);
        assert(m>0);
        break;

    }

    src = (attr_id_t *) malloc (m * sizeof(attr_id_t));
    dest = (attr_id_t *) malloc(m * sizeof(attr_id_t));
    degree = (attr_id_t *) calloc(n, sizeof(attr_id_t));

    assert(src != NULL);
    assert(dest != NULL);
    assert(degree != NULL);

    int_weight = (int *) malloc(m * sizeof(int));
    assert(int_weight != NULL);

    count = 0;

    while (fgets(buf, 500, fp) != NULL)  {
        sscanf(buf, "%*c %ld %ld %d", &u, &v, &int_wt);
        if ((u <= 0) || (u > n) || (v <= 0) || (v > n)) {
            fprintf(stderr, "Error reading edge # %d in the input file."
                    " Please check the input graph file.\n", count+1);
            exit(1);
        }
        src[count] = u-1;
        dest[count] = v-1;
        degree[u-1]++;
        degree[v-1]++;
        int_weight[count] = int_wt;
        count++;
    }

    fclose(fp);

    if (count != m) {
        fprintf(stderr, "Error! Number of edges specified in problem line (%ld)" 
                " does not match the total number of edges (%ld) in file."
                " Please check"
                " the graph input file.\n", m, count);
        exit(1);
    }

    /*
       for (i=0; i<m; i++) {
       fprintf(stderr, "[%d %d] ", src[i], dest[i]);
       }
     */

    G->endV = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
    assert(G->endV != NULL);

    G->edge_id = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
    assert(G->edge_id != NULL);

    G->numEdges = (attr_id_t *) malloc((n+1)*sizeof(attr_id_t));
    assert(G->numEdges != NULL);

    G->undirected = 1;
    G->weight_type = 1;
    G->zero_indexed = 0;

    G->n = n;
    G->m = 2*m;

    G->int_weight_e = (int *) malloc(G->m * sizeof(int));       
    assert(G->int_weight_e != NULL);

    /* ToDo: parallelize this step */
    G->numEdges[0] = 0; 
    for (i=1; i<=G->n; i++) {
        G->numEdges[i] = G->numEdges[i-1] + degree[i-1];
    }

    for (i=0; i<count; i++) {
        u = src[i];
        v = dest[i];
        offset = degree[u]--;
        G->endV[G->numEdges[u]+offset-1] = v;
        G->int_weight_e[G->numEdges[u]+offset-1] = int_weight[i];
        G->edge_id[G->numEdges[u]+offset-1] = i;

        offset = degree[v]--;
        G->endV[G->numEdges[v]+offset-1] = u;
        G->int_weight_e[G->numEdges[v]+offset-1] = int_weight[i];
        G->edge_id[G->numEdges[v]+offset-1] = i;
    } 

    /*          
    for (i=0; i<G->n; i++) {
        for (j=G->numEdges[i]; j<G->numEdges[i+1]; j++) {
            fprintf(stderr, "<%ld %ld %d> ", i+1, G->endV[j]+1, 
            G->int_weight_e[j]);
        }
    }
    */

    free(buf);
    free(degree);
    free(src);
    free(dest);
}
