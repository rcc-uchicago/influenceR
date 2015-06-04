#include "graph_defs.h"
#include "graph_gen.h"

void read_SNAP_graph(graph_t* G, char* filename) {

    char *buf;
    FILE *fp;

    long i;
    long undirected, repeated, zero_indexed, weight_type;
    char udc, zic, wtc;
    long n, m; 
    long count, offset;
    int int_wt, *int_weight;
    long l_wt, *l_weight;
    float fl_wt, *fl_weight;
    double dbl_wt, *dbl_weight;
    long u, v;
    attr_id_t *src;
    attr_id_t *dest;
    attr_id_t *degree;
    count=0;

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

    while (fgets(buf, 500, fp) != NULL)  {

        /* Skip all lines until we reach the problem line */
        if (buf[0] != 'p')
            continue;

        sscanf(buf, "%*c %ld %ld %c %c %c", &n, &m, &udc, &wtc,
                &zic);

        assert(n>0);
        assert(m>0);
        assert((udc == 'u') || (udc == 'd') || (udc == 'r'));
        assert((wtc == 'u') || (wtc == 'i') || (wtc == 'l') || (wtc == 'f') ||
                (wtc == 'd'));
        assert((zic == '0') || (zic == '1'));

        repeated = 0;
        if (udc == 'u') {
            undirected = 1;
        } else if (udc == 'd') {
            undirected = 0;
        } else if (udc == 'r') {
            undirected = 1;
            repeated   = 1;
        }

        if (wtc == 'u') {
            weight_type = 0;
        } else if (wtc == 'i') {
            weight_type = 1;
        } else if (wtc == 'l') {
            weight_type = 2;
        } else if (wtc == 'f') {
            weight_type = 3;
        } else {
            weight_type = 4;
        }

        if (zic == '0') {
            zero_indexed = 1;
        } else {
            zero_indexed = 0;
        }

        break;
    }

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

    count = 0;

    if ((zero_indexed == 1) && (undirected == 0)) {
        if (weight_type == 0) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld", &u, &v);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                count++;
            }
        }

        if (weight_type == 1) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %d", &u, &v, &int_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                int_weight[count] = int_wt;
                count++;
            }
        }

        if (weight_type == 2) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %ld", &u, &v, &l_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                l_weight[count] = l_wt;
                count++;
            }
        }

        if (weight_type == 3) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %f", &u, &v, &fl_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                fl_weight[count] = fl_wt;
                count++;
            }
        }

        if (weight_type == 4) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %lf", &u, &v, &dbl_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                dbl_weight[count] = dbl_wt;
                count++;
            }
        }
    }


    if ((zero_indexed == 0) && (undirected == 0)) {
        if (weight_type == 0) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld", &u, &v);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                count++;
            }
        }

        if (weight_type == 1) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %d", &u, &v, &int_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                int_weight[count] = int_wt;
                count++;
            }
        }

        if (weight_type == 2) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %ld", &u, &v, &l_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                l_weight[count] = l_wt;
                count++;
            }
        }

        if (weight_type == 3) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %f", &u, &v, &fl_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                fl_weight[count] = fl_wt;
                count++;
            }
        }

        if (weight_type == 4) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %lf", &u, &v, &dbl_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                dbl_weight[count] = dbl_wt;
                count++;
            }
        }
    }


    if ((zero_indexed == 1) && (undirected == 1)) {
        if (weight_type == 0) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld", &u, &v);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                if ((repeated) && (u > v)) continue;
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                degree[v]++;
                count++;
            }
        }

        if (weight_type == 1) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %d", &u, &v, &int_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                if ((repeated) && (u > v)) continue;
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                degree[v]++;
                int_weight[count] = int_wt;
                count++;
            }
        }

        if (weight_type == 2) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %ld", &u, &v, &l_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                if ((repeated) && (u > v)) continue;
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                degree[v]++;
                l_weight[count] = l_wt;
                count++;
            }
        }

        if (weight_type == 3) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %f", &u, &v, &fl_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                if ((repeated) && (u > v)) continue;
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                degree[v]++;
                fl_weight[count] = fl_wt;
                count++;
            }
        }

        if (weight_type == 4) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %lf", &u, &v, &dbl_wt);
                assert((u>=0)&&(u<n)&&(v>=0)&&(v<n));
                if ((repeated) && (u > v)) continue;
                src[count] = u;
                dest[count] = v;
                degree[u]++;
                degree[v]++;
                dbl_weight[count] = dbl_wt;
                count++;
            }
        }
    }


    if ((zero_indexed == 0) && (undirected == 1)) {
        if (weight_type == 0) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld", &u, &v);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                if ((repeated) && (u > v)) continue;
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                degree[v-1]++;
                count++;
            }
        }

        if (weight_type == 1) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %d", &u, &v, &int_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                if ((repeated) && (u > v)) continue;
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                degree[v-1]++;
                int_weight[count] = int_wt;
                count++;
            }
        }

        if (weight_type == 2) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %ld", &u, &v, &l_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                if ((repeated) && (u > v)) continue;
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                degree[v-1]++;
                l_weight[count] = l_wt;
                count++;
            }
        }

        if (weight_type == 3) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %f", &u, &v, &fl_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                if ((repeated) && (u > v)) continue;
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                degree[v-1]++;
                fl_weight[count] = fl_wt;
                count++;
            }
        }

        if (weight_type == 4) {
            while (fgets(buf, 500, fp) != NULL)  {
                sscanf(buf, "%ld %ld %lf", &u, &v, &dbl_wt);
                assert((u>0)&&(u<=n)&&(v>0)&&(v<=n));
                if ((repeated) && (u > v)) continue;
                src[count] = u-1;
                dest[count] = v-1;
                degree[u-1]++;
                degree[v-1]++;
                dbl_weight[count] = dbl_wt;
                count++;
            }
        }
    }

    fclose(fp);

    if ((repeated == 0) && (count != m)) {
        fprintf(stderr, "Error! Number of edges specified in problem line (%ld)" 
                " does not match the total number of edges (%ld) in file." 
                " Please check"
                " the graph input file. Exiting ...\n", m, count);
        exit(1);
    }

    /*
    for (i=0; i<m; i++) {
        fprintf(stderr, "[%d %d] ", src[i], dest[i]);
    }
     */
    if (repeated) m = count;

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

    G->undirected = undirected;
    G->weight_type = weight_type;
    G->zero_indexed = zero_indexed;

    G->n = n;
    if (G->undirected) {
        G->m = 2*m;
    }
    else {
        G->m = m;
    }

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

    /* ToDo: parallelize this step */
    G->numEdges[0] = 0; 
    for (i=1;i<=G->n;i++) {
        G->numEdges[i] = G->numEdges[i-1] + degree[i-1];
    }

    for (i=0; i<count; i++) {
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

    /*       
    for (i=0; i<G->n; i++) {
        for (j=G->numEdges[i]; j<G->numEdges[i+1]; j++) {
            fprintf(stderr, "<%ld %ld> ", i, G->endV[j]);
        }
    }
    */

    if (weight_type == 1) {
	free(int_weight);
    } else if (weight_type == 2) {
	free(l_weight);
    } else if (weight_type == 3) {
	free(fl_weight);
    } else if (weight_type == 4) {
	free(dbl_weight);
    } 

    free(buf);
    free(degree);
    free(src);
    free(dest);
}
