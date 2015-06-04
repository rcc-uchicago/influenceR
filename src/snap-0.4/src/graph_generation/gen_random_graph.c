#include "graph_defs.h"
#include "graph_gen.h"
#include "sprng.h"

void gen_random_graph(graph_t* G, char* filename) {

    long i;
    long undirected, weight_type;
    long n, m; 
    long offset;
    int *int_weight;
    long *l_weight;
    float *fl_weight;
    double *dbl_weight;
    double min_weight, max_weight;
    long u, v;
    int* stream, seed;
    attr_id_t *src;
    attr_id_t *dest;
    attr_id_t *degree;

    read_random_config_file(G, filename);

    undirected = G->undirected;
    weight_type = G->weight_type;
    n = G->n;
    m = G->m;

    src = (attr_id_t *) malloc (m * sizeof(attr_id_t));
    dest = (attr_id_t *) malloc(m * sizeof(attr_id_t));
    degree = (attr_id_t *) calloc(n, sizeof(attr_id_t));

    assert(src != NULL);
    assert(dest != NULL);
    assert(degree != NULL);

    if (weight_type == 1) {
        int_weight = (int *) malloc(m * sizeof(int));
        assert(int_weight != NULL);
    } else if (weight_type == 2) {
        l_weight = (long *) malloc(m * sizeof(long));
        assert(l_weight != NULL);
    } else if (weight_type == 3) {
        fl_weight = (float *) malloc(m * sizeof(float));
        assert(fl_weight != NULL);
    } else if (weight_type == 4) {
        dbl_weight = (double *) malloc(m * sizeof(double));
        assert(dbl_weight != NULL);
    }

    /* Initialize RNG stream */
    seed = 2387;
    stream = init_sprng(0, 0, 1, seed, SPRNG_DEFAULT);

    /* Generate edges */
    for (i=0; i<m; i++) {
        u = n * sprng(stream);
        v = n * sprng(stream);
        src[i] = u;
        dest[i] = v;
        degree[u]++;
        if (undirected)
            degree[v]++;
    }

    min_weight = G->min_weight;
    max_weight = G->max_weight;

    /* Generate edge weights */
    if (weight_type == 1) {
        for (i=0; i<m; i++) {
            int_weight[i]  = min_weight + 
                (int) ((max_weight-min_weight)*sprng(stream));
        }
    } else if (weight_type == 2) {
        for (i=0; i<m; i++) {
            l_weight[i]  = min_weight + 
                (long) (max_weight-min_weight)*sprng(stream);
        }
    } else if (weight_type == 3) {
        for (i=0; i<m; i++) {
            fl_weight[i]  = min_weight + 
                (float) (max_weight-min_weight)*sprng(stream);
        }
    } else if (weight_type == 4) {
        for (i=0; i<m; i++) {
            dbl_weight[i]  = min_weight + 
                (max_weight-min_weight)*sprng(stream);
        }
    } 

    /* Update graph data structure */
    if (undirected) {
        G->endV = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
        G->edge_id = (attr_id_t *) malloc(2*m * sizeof(attr_id_t));
    } else {
        G->endV = (attr_id_t *) calloc(m, sizeof(attr_id_t));
        G->edge_id = (attr_id_t *) malloc(m * sizeof(attr_id_t));   
    }
    assert(G->endV != NULL);
    assert(G->edge_id != NULL);

    G->numEdges = (attr_id_t *) malloc((n+1)*sizeof(attr_id_t));
    assert(G->numEdges != NULL);

    G->n = n;
    if (G->undirected)
        G->m = 2*m;
    else
        G->m = m;

    if (weight_type) {

        if (weight_type == 1) {
            G->int_weight_e = (int *) malloc(G->m * sizeof(int));
            assert(G->int_weight_e != NULL);
        }

        if (weight_type == 2) {
            G->l_weight_e = (long *) malloc(G->m * sizeof(long));
            assert(G->l_weight_e != NULL);
        }

        if (weight_type == 3) {
            G->fl_weight_e = (float *) malloc(G->m * sizeof(float));
            assert(G->fl_weight_e != NULL);
        }

        if (weight_type == 4) {
            G->dbl_weight_e = (double *) malloc(G->m * sizeof(double));       
            assert(G->dbl_weight_e != NULL);
        }

    }

    G->numEdges[0] = 0; 
    for (i=1;i<=G->n;i++) {
        G->numEdges[i] = G->numEdges[i-1] + degree[i-1];
    }

    for (i=0; i<m; i++) {
        u = src[i];
        v = dest[i];
        offset = degree[u]--;
        G->endV[G->numEdges[u]+offset-1] = v;
        G->edge_id[G->numEdges[u]+offset-1] = i;
        if (weight_type) {
            if (weight_type == 1) {
                G->int_weight_e[G->numEdges[u]+offset-1] = int_weight[i];
            }
            if (weight_type == 2) {
                G->l_weight_e[G->numEdges[u]+offset-1] = l_weight[i];
            }
            if (weight_type == 3) {
                G->fl_weight_e[G->numEdges[u]+offset-1] = fl_weight[i];
            }
            if (weight_type == 4) {
                G->dbl_weight_e[G->numEdges[u]+offset-1] = dbl_weight[i];
            }
        }

        if (undirected) {
            offset = degree[v]--;
            G->endV[G->numEdges[v]+offset-1] = u;
            G->edge_id[G->numEdges[v]+offset-1] = i;
            if (weight_type) {
                if (weight_type == 1) {
                    G->int_weight_e[G->numEdges[v]+offset-1] = int_weight[i];
                }
                if (weight_type == 2) {
                    G->l_weight_e[G->numEdges[v]+offset-1] = l_weight[i];
                }
                if (weight_type == 3) {
                    G->fl_weight_e[G->numEdges[v]+offset-1] = fl_weight[i];
                }
                if (weight_type == 4) {
                    G->dbl_weight_e[G->numEdges[v]+offset-1] = dbl_weight[i];
                }
            }

        }
    }
    G->zero_indexed = 1;
    
    /*
    for (i=0; i<G->n; i++) {
        for (j=G->numEdges[i]; j<G->numEdges[i+1]; j++) {
            fprintf(stderr, "<%ld %ld> ", i, G->endV[j]);
        }
    }
    */

    free(src);
    free(dest);
    free(degree);

    if (weight_type == 1) 
        free(int_weight);
    if (weight_type == 2)
        free(l_weight);
    if (weight_type == 3)
        free(fl_weight);
    if (weight_type == 4)
        free(dbl_weight);

    free_sprng(stream);
}

void read_random_config_file(graph_t* G, char* configfile) {

    /* read parameters from config file */
    FILE *fp;
    char line[128], var[32];
    double val;

    fp = fopen(configfile,"r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open config file: %s\n",configfile);
        exit(-1);
    }

    while (fgets(line, sizeof (line), fp) != NULL) {
        sscanf(line, "%s %lf", var, &val);
        if (*var == '#') continue;  /* comment */
        if (strcmp(var, "n") == 0) {
            G->n = (long) val;
            assert(G->n > 0);
        } else if (strcmp(var, "m") == 0) {
            G->m = (long) val;
            assert(G->m > 0);
        } else if (strcmp(var, "undirected") == 0) {
            G->undirected = (int) val;
            assert((G->undirected == 0) || (G->undirected == 1)); 
        } else if (strcmp(var, "weight_type") == 0) {
            G->weight_type = (int) val;
            assert(G->weight_type >= 0);
            assert(G->weight_type <= 4);
        } else if (strcmp(var, "max_weight") == 0) {
            G->max_weight = val;
        } else if (strcmp(var, "min_weight") == 0) {
            G->min_weight = val;
        } else {
            fprintf(stderr,"Unknown parameter: %s\n", line);
        }
    }

    fclose(fp);

}

