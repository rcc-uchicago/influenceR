#include "graph_defs.h"

void adj_bheap_insert(adj_bheap_t *h, attr_id_t val)
{

    attr_id_t i, j, y;

    i = ++(h->count);
    while (i >= 2) {
        j = i / 2;
        y = h->vals[j];
        if (val >= y) 
            break;
        h->vals[i] = y;
        i = j;
    }
    h->vals[i] = val;
}


