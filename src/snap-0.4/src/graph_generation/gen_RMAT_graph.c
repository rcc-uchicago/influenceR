#include <time.h>
#include <sys/types.h>
#include <unistd.h>

#include "graph_defs.h"
#include "graph_gen.h"
#include "sprng.h"

void gen_RMAT_graph(graph_t* G, char* filename) {

    long i, j;
    long undirected, weight_type;
    long n, m; 
    long offset;
    double *params, a, b, c, d;
    double av, bv, cv, dv, S, p;
    int SCALE;
    double var;
    long step;
    int *int_weight;
    long *l_weight;
    float *fl_weight;
    double *dbl_weight;
    double min_weight, max_weight;
    long permute_vertices;
    attr_id_t* permV, tmpVal;
    long u, v;
    int* stream, seed;
    attr_id_t *src;
    attr_id_t *dest;
    attr_id_t *degree;

    params = (double *) malloc(4*sizeof(double));
    read_RMAT_config_file(G, filename, params);
    a = params[0];
    b = params[1];
    c = params[2];
    assert(a+b+c < 1);
    d = 1  - (a+b+c);

    permute_vertices = (long) params[3];

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
    if (getenv ("SEED_WITH_TIME"))
      seed = ((int) time (NULL)) ^ ((int) getpid ());
    else
      seed = 2387;
    stream = init_sprng(0, 0, 1, seed, SPRNG_DEFAULT);
    SCALE = log2(n);
    fprintf(stderr, "Scale: %d\n", SCALE);

    /* Generate edges */
    for (i=0; i<m; i++) {

        u = 1;
        v = 1;
        step = n/2;

        av = a;
        bv = b;
        cv = c;
        dv = d;

        p = sprng(stream);
        if (p < av) {
            /* Do nothing */
        } else if ((p >= av) && (p < av+bv)) {
            v += step;
        } else if ((p >= av+bv) && (p < av+bv+cv)) {
            u += step;
        } else {
            u += step;
            v += step;
        }

        for (j=1; j<SCALE; j++) {
            step = step/2;

            /* Vary a,b,c,d by up to 10% */
            var = 0.1;
            av *= 0.95 + var * sprng(stream);
            bv *= 0.95 + var * sprng(stream);
            cv *= 0.95 + var * sprng(stream);
            dv *= 0.95 + var * sprng(stream);

            S = av + bv + cv + dv;
            av = av/S;
            bv = bv/S;
            cv = cv/S;
            dv = dv/S;

            /* Choose partition */
            p = sprng(stream);
            if (p < av) {
                /* Do nothing */
            } else if ((p >= av) && (p < av+bv)) {
                v += step;
            } else if ((p >= av+bv) && (p < av+bv+cv)) {
                u += step;
            } else {
                u += step;
                v += step;
            }
        }

        src[i] = u-1;
        dest[i] = v-1;
    }

    if (permute_vertices) {

        permV = (attr_id_t *) malloc(n*sizeof(attr_id_t));
        assert(permV != NULL);

        for (i=0; i<n; i++) {
            permV[i] = i;
        }

        for (i=0; i<n; i++) {
            j = n * sprng(stream);
            tmpVal = permV[i];
            permV[i] = permV[j];
            permV[j] = tmpVal;
        }

        for (i=0; i<m; i++) {
            src[i]  = permV[src[i]];
            dest[i] = permV[dest[i]];
        }

    }

    for (i=0; i<m; i++) {
        degree[src[i]]++;
        if (undirected)
            degree[dest[i]]++;
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
                (float) floor((max_weight-min_weight)*sprng(stream));
        }
    } else if (weight_type == 4) {
        for (i=0; i<m; i++) {
            dbl_weight[i]  = min_weight + 
                floor((max_weight-min_weight)*sprng(stream));
        }
    } 

    /* Update graph data structure */
    if (undirected) {
        G->endV = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
        G->edge_id = (attr_id_t *) calloc(2*m, sizeof(attr_id_t));
    } else { 
        G->endV = (attr_id_t *) calloc(m, sizeof(attr_id_t));
        G->edge_id = (attr_id_t *) calloc(m, sizeof(attr_id_t));
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
    free(params);

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

void read_RMAT_config_file(graph_t* G, char* configfile, double* params) {

    /* read parameters from config file */
    FILE *fp;
    char line[128], var[32];
    double val;

    fp = fopen(configfile,"r");
    if (fp == NULL) {
        fprintf(stderr, "Unable to open config file:i %s\n",configfile);
        exit(-1);
    }

    /* default values */
    params[0] = 0.45; params[1] = 0.25; params[2] = 0.15; params[3] = 1; 
    while (fgets(line, sizeof (line), fp) != NULL) {
        sscanf(line, "%s %lf", var, &val);
        if (*var == '#') continue;  /* comment */
        if (strcmp(var, "n") == 0) {
            G->n = (long) val;
            assert(G->n > 0);
        } else if (strcmp(var, "m") == 0) {
            G->m = (long) val;
            assert(G->m > 0);
        } else if (strcmp(var, "a") == 0) {
            params[0] = val;
            assert((val > 0) && (val < 1));
        } else if (strcmp(var, "b") == 0) {
            params[1] = val;
            assert((val > 0) && (val < 1));
        } else if (strcmp(var, "c") == 0) {
            params[2] = val;
            assert((val > 0) && (val < 1)); 
        } else if (strcmp(var, "permute_vertices") == 0) {
            params[3] = val;
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

