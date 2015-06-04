#include "graph_defs.h"
#include "graph_kernels.h"

struct s_ent {
    attr_id_t a, b;
};

struct s_ent *BiCC_stack;

attr_id_t *dfsnumber, *highwater, lastdfsnumber=0, top=-1, count=0, *cc;
attr_id_t *color, *Low, *d, *pred, *art, art_count;
attr_id_t timex;

attr_id_t biconnected_components(graph_t* G, attr_id_t* component_num) {

    long i;
    long n, m;
    attr_id_t bcc_total;

    n = G->n;
    m = G->m;

    dfsnumber = (attr_id_t *) malloc(n*sizeof(attr_id_t));
    assert(dfsnumber != NULL);
    highwater = (attr_id_t *) malloc(n*sizeof(attr_id_t));
    assert(highwater != NULL);
#if 0
    cc = (attr_id_t *) malloc(n*sizeof(attr_id_t));
#endif
    BiCC_stack = (struct s_ent *) malloc(m*sizeof(struct s_ent));
    assert(BiCC_stack != NULL);

    cc = component_num;
    for (i=0; i<n; i++) {
        dfsnumber[i]=0;
        cc[i]=0;
    }

    biconnected_components_recursive(G, 0, -1);

    free(dfsnumber);
    free(highwater);
    free(BiCC_stack);
#if 0
    free(cc);
#endif

    bcc_total = -1;
    cc = component_num;
    for (i=0 ; i<n ; i++)
        if (cc[i] > bcc_total)
            bcc_total = cc[i];
    bcc_total++;
    assert(bcc_total > 0);

    return(bcc_total);
}

void find_articulation_points(graph_t* G, attr_id_t* art_point_map) {

    long n;
    n = G->n;

    color = (attr_id_t *) calloc(n, sizeof(attr_id_t));
    art   = (attr_id_t *) calloc(n, sizeof(attr_id_t));
    Low   = (attr_id_t *) malloc(n*sizeof(attr_id_t));
    pred  = (attr_id_t *) malloc(n*sizeof(attr_id_t)); 
    d     = (attr_id_t *) malloc(n*sizeof(attr_id_t));

    timex = 0;
    pred[0] = -1;
    art_count = 0; 

    art_points_recursive(G, 0);

    free(color);
    free(art);
    free(Low);
    free(pred);
    free(d);

}

void art_points_recursive(graph_t* G, attr_id_t u) {
    long i;
    attr_id_t v;

    color[u] = 1;
    Low[u] = d[u] = ++timex;

    for (i=G->numEdges[u]; i<G->numEdges[u+1]; i++) {
        v = G->endV[i];
        if (color[v] == 0) {
            pred[v] = u;
            art_points_recursive(G, v);
            if (Low[v] < Low[u]) {
                Low[u] = Low[v];
            }
            if (pred[u] == -1) {
                if ((G->numEdges[u+1] - G->numEdges[u]) > 1) {
                    art[u] = 1;
                }
            } else if (Low[v] >= d[u]) {
                art[u] = 1;
            }
        } else if (v != pred[u]) {
            if (d[v] < Low[u]) {
                Low[u] = d[v];
            }
        }
    } 

}

void biconnected_components_recursive(graph_t* G, attr_id_t v, attr_id_t u) {

    attr_id_t i,w,a,b;

    lastdfsnumber++;
    dfsnumber[v]=lastdfsnumber;
    highwater[v]=lastdfsnumber;

    for (i=G->numEdges[v]; i<G->numEdges[v+1]; i++) {
        w = G->endV[i];
        if (dfsnumber[w]==0) {
            BiCC_stack_push(v,w);
            biconnected_components_recursive(G, w, v);
            if (highwater[w] < highwater[v]) 
                highwater[v] = highwater[w];
            if (dfsnumber[v] <= highwater[w]) {
                count++;
                BiCC_stack_pop(&a,&b);
                cc[a]=count;
                cc[b]=count;
                while (!(a==v && b==w)) {
                    BiCC_stack_pop(&a, &b);
                    if (cc[a] == 0) 
                        cc[a] = count;
                    if (cc[b] == 0) 
                        cc[b] = count;
                }
            } 
        } else if((dfsnumber[w] < dfsnumber[v]) && (w != u)) {
            BiCC_stack_push(v,w);
            if (dfsnumber[w] < highwater[v]) 
                highwater[v] = dfsnumber[w];
        }
    }

}


void BiCC_stack_push(attr_id_t a, attr_id_t b) {
    top++;
    BiCC_stack[top].a=a;
    BiCC_stack[top].b=b;
}

void BiCC_stack_pop(attr_id_t* a, attr_id_t* b) {
    *a = BiCC_stack[top].a;
    *b = BiCC_stack[top].b;
    top--;
}
