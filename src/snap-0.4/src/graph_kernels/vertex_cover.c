#include "graph_kernels.h"

int vertex_cover_weighted(graph_t *G)
{    

    double *wp_v;                /* weight associated with each vertex. */
    double *delta_e;             /* Delta of each edge as mentioned in 
                                    the algorithm */
    attr_id_t *degree_v;                /* Degree of each vertex */
    attr_id_t *visited_e, *visited_v;   /* Whether that edge has been 
                    visited or not. Ditto for vertex. Visited_v is the 
                    final cover. */
    attr_id_t *position_e;       /* Position stores the corresponding 
                                    position of the undirected edge for 
                                    each edge.
                                    */
    attr_id_t i,j,u,v,n,k, edge_counter;

    double *memblock;
    attr_id_t *memblock1;
    double val1,val2;
    int count;
    double sum;

    memblock = (double*) malloc(sizeof(double)*(G->n+2*G->m));
    wp_v = memblock;
    delta_e = memblock + G->n;
    memblock1 = (attr_id_t*) malloc(sizeof(attr_id_t)*(2*G->n + 4*G->m));
    degree_v = memblock1;
    visited_v = memblock1 + G->n;
    visited_e = memblock1 + 2*G->n;
    position_e = memblock1 + 2*(G->n + G->m);
    n = G->n;
#ifdef _OPENMP
    #pragma omp parallel for private(u,v,j,k)
#endif
    for(i=0; i<n; i++)
    {
        wp_v[i] = G->dbl_weight_v[i];
        degree_v[i] = G->numEdges[i+1] - G->numEdges[i];
        visited_v[i] = 0;
        if(degree_v[i] == 0)
            visited_v[i]=1;
        for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
        {
            u = i;
            v = G->endV[j];
            delta_e[j] = 0;
            visited_e[j] = 0;
            if(v < u )
                continue;        /* we have already covered this case 
                                    when we visited v. */
            for (k=G->numEdges[v]; k<G->numEdges[v+1]; k++)
            {
                if(G->endV[k] == u)
                    break;
            }
            position_e[j] = k;
            position_e[k] = j;
        }
    }

    edge_counter = 2*G->m;
    count =0;
    while(edge_counter > 0)
    {
        count ++;
#ifdef _OPENMP
        #pragma omp parallel
        {
        #pragma omp for private(j,u,v,val1,val2)
#endif
        for (i=0; i< n; i++)
        {
            if (visited_v[i] == 1)
                continue;
            for(j=G->numEdges[i]; j< G->numEdges[i+1]; j++)
            {
                if(visited_e[j] == 1 )
                    continue;
                u = i;
                v = G->endV[j];
                val1 = wp_v[u]/degree_v[u];
                val2 = wp_v[v]/degree_v[v];
                delta_e[j] = val1 < val2 ? val1 : val2;
            }
        }
#ifdef _OPENMP
        #pragma omp for private(j) reduction(-:edge_counter)
#endif
        for(i=0; i<n; i++)
        {
            if (visited_v[i] == 1)
                continue;
            sum = 0.0;
            for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
            {
                if(visited_e[j] == 1)
                    continue;
                sum += delta_e[j];
            }
            wp_v[i] -= sum;
            if(wp_v[i] <= 0.00001)    /* aka this vertex is in VC. */
            {
                visited_v[i] = 1;
                edge_counter -= degree_v[i]*2;  /* It is multiplied by 
                                                   because it is an 
                                                   undirected graph.
                                                 */
                for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
                {
                    if(visited_e[j] == 1)
                        continue;
                    visited_e[j] = 1;
                    degree_v[G->endV[j]] -= 1;
                    visited_e[position_e[j]] = 1;
                }
                
            }
        }
        }

#ifdef _OPENMP        
    }
#endif

    count=0;
    for(i=0; i<G->n; i++)
    {
        if(visited_v[i] == 1 && ((G->numEdges[i+1] - G->numEdges[i])>0))
            count ++;
    }
    free(memblock);
    free(memblock1);
    return count;
}




int vertex_cover_unweighted(graph_t *G)
{
    attr_id_t i,j,u,v, n, m;
    attr_id_t max, max_e, max_u, max_v,edge_counter;
    
    attr_id_t *visited_v, *visited_e;
    attr_id_t * degree_v;
    attr_id_t count;
    attr_id_t *memblock;


    memblock = (attr_id_t*) malloc(sizeof(attr_id_t)*(2*G->m + 2*G->n));
    visited_v = memblock;
    degree_v = memblock + G->n;
    visited_e = memblock + 2*G->n;
    n = G->n;
    m = G->m;
#ifdef _OPENMP    
    #pragma omp parallel 
    {
        #pragma omp for
#endif
        for(i=0; i<n; i++)
        {
            visited_v[i] = 0;
            degree_v[i] = G->numEdges[i+1] - G->numEdges[i];
        }
#ifdef _OPENMP
        #pragma omp for
#endif
        for(i=0; i<2*m; i++)
        {
            visited_e[i] = 0;
        }
#ifdef _OPENMP
    }
#endif
    edge_counter = 2*G->m;
    while(edge_counter > 0)
    {
        max = 0;
#ifdef _OPENMP
        #pragma omp parallel for shared(max,max_e,max_u,max_v) private(j,u,v)
#endif
        for(i=0; i<n; i++)
        {
            if(degree_v[i] == 0)
                continue;
            for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
            {
                u = i;
                v = G->endV[j];
                if(degree_v[u] + degree_v[v] > max)
                {
                    max = degree_v[u]+ degree_v[v];
                    max_e = j;
                    max_u = u;
                    max_v = v;
                }
            }
        }
        edge_counter -= max;
        visited_e[max_e] =1;
        degree_v[max_u] = 0;
        visited_v[max_u] = 1;
        degree_v[max_v] = 0;
        visited_v[max_v] = 1;

    }
    count = 0;
    for(i=0; i<G->n; i++)
    {
        if(visited_v[i] == 1)
            count ++;
    }
    free(memblock);
    return count;

}
