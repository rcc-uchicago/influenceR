#include <time.h>

#include "graph_partitioning.h"
#include "utils.h"

typedef struct {
    double *dist;   /* stores the value */
    int *vertex;    /* map from the index of dist to vertex */
    int *map;       /* map from vertex to index of dist */
    int *from;      /* map from vertex to its parents in djkstra */
    int size;
} heap;

void computeEigen(graph_t *G, double *eigenVectorOld, 
        double *eigenVectorNew, attr_id_t *v2C, attr_id_t *v2pos, 
        attr_id_t* degree, attr_id_t *vertex, attr_id_t *toSplit, 
        attr_id_t currCommunity, attr_id_t communitySize, 
        attr_id_t degreeSum);

void computeModularityValue(graph_t *G, attr_id_t *membership, 
        attr_id_t numCommunities, double *modularity) {
    attr_id_t i,j;
    attr_id_t n, m;
    attr_id_t comm;
    double mod=0.0;
    double degree_u,degree_v;
    n = G->n;
    m = G->m; 

    for(i=0; i<n; i++)
    {
        comm = membership[i];
        degree_u = (double)G->numEdges[i+1] - G->numEdges[i];
        for(j=G->numEdges[i]; j<G->numEdges[i+1]; j++)
        {
            if(membership[G->endV[j]] == comm)
                mod +=1.0;
        }
        for(j=0; j<G->n; j++)
        {
            degree_v = (double)G->numEdges[j+1] - G->numEdges[j];
            if(comm == membership[j])
                mod -= (double)(degree_u*degree_v)/(double)(2.0*G->m);
        }
    }
    *modularity = mod/(2*G->m);
}


void modularity_spectral(graph_t *G, attr_id_t *membership, 
        attr_id_t *numCommunities, attr_id_t use_improvement) {
    attr_id_t *v2C, *degree, *vertex, *v2pos;
    attr_id_t u,v,curCommunity=0, newCommunity,toSplit;
    attr_id_t sumV1,sumV2,comm,count1,count2;
    attr_id_t i,j,communitySize,degreeSum;
    list_t *Q;
    node_t *first;
    double *eigenVectorOld,*eigenVectorNew;
    attr_id_t continue_flag = 0;    
    double *contribution, max_contrib, modularity, new_modularity;
    attr_id_t degree_u, degree_v, degreeSum1, degreeSum2, maxv, flag, flag_counter;
    attr_id_t n, m;

    n = G->n;
    m = G->m;
    *numCommunities = 1;

    v2C = membership;         /* v2C is a map from vertex to the 
                                 community it belongs. */
    vertex = (attr_id_t *) malloc(sizeof(attr_id_t)*n);  
    /* vertex is an array of vertices belonging to a particular community.*/

    v2pos  = (attr_id_t *) malloc(sizeof(attr_id_t)*n);  
    /* v2pos is a map from vertex to its respective position in the community.*/

    degree =  (attr_id_t *) malloc(sizeof(attr_id_t)*n); 
    /* degree is a map from vertex to its degree in its respective community.*/

    eigenVectorOld=(double*)malloc(sizeof(double)*n);
    eigenVectorNew=(double*)malloc(sizeof(double)*n);
    contribution = (double*)malloc(sizeof(double)*n);

    assert(v2C != NULL); assert(vertex != NULL); assert(degree != NULL);
    assert(eigenVectorOld!=NULL); assert(eigenVectorNew!=NULL);
    assert(contribution != NULL);

    for(i=0; i< G->n; i++)
    {
        v2C[i] = 0;
        vertex[i] = i;
        v2pos[i] = i;
        degree[i] = G->numEdges[i+1] - G->numEdges[i];
    }
    /* Making a queue. This queue will store all the 
       communities that are yet to be processed. */
    Q=(list_t*)makeList();
    append(Q, makeNode(curCommunity));
    /* printList(Q); */


    while(Q->size > 0)
    {
        first = (node_t*) getFirst(Q);
        curCommunity = first->id;
        deleteFirst(Q);
        continue_flag = 0;
#ifdef _OPENMP
OMP("omp parallel for")
#endif
        for(i=0; i<n; i++)
        {
            vertex[i] = -1;
        }
        degreeSum=0;
        communitySize=-1;
        /* Checking which all vertices belong to this community 
           and updating the vertex Vector accordingly. */
        for(i=0; i<G->n; i++)
        {
            contribution[i] = 0.0;        /* added later for klin */
            if(v2C[i] == curCommunity)
            {
                {
                    communitySize++;
                    vertex[communitySize] = i;
                }
                v2pos[i] = communitySize;
                degreeSum += G->numEdges[i+1]-G->numEdges[i];
            }
        }

        communitySize ++;
        if(communitySize == 1)    continue;

        /* Calculating modularity by the current Community. */   
        modularity = 0.0;
#ifdef _OPENMP
OMP("omp parallel for private(j) reduction(+:modularity)")
#endif
        for(i=0; i<n;i++)
        {
            if(v2C[i] == curCommunity)
            {
                for(j=G->numEdges[i] ; j<G->numEdges[i+1]; j++)
                {
                    if(v2C[G->endV[j]] == curCommunity)
                        modularity += 1.0;
                }
                modularity += 
                    -((G->numEdges[i+1] - G->numEdges[i])*degreeSum)/(2.0*G->m);
            }
        }
        modularity /=(2.0*G->m);

        /* Computing eigen vector. */
        computeEigen(G,eigenVectorOld, eigenVectorNew, v2C, v2pos, degree, 
                vertex, &toSplit, curCommunity, communitySize, degreeSum);
        if(toSplit == 0)
            continue;

        newCommunity = *numCommunities;
        count1=count2=sumV1=sumV2=0;
        degreeSum1 = degreeSum2 = 0;
        new_modularity=0.0;
#ifdef _OPENMP    
OMP("omp parallel for reduction(+:count1,count2)")
#endif
        for(i=0; i<communitySize; i++)
        {
            if(eigenVectorOld[i] > 0) count1++;
            else count2++;
        }
        if(count1 == 0 || count2 == 0)
            continue;          /* All eigen values are of same size 
                                  and hence no division is required. */

        /* Now, we actually divide the community to new communities. */
#ifdef _OPENMP
OMP("omp parallel if (communitySize>100) ")
#endif
        {
#ifdef _OPENMP
OMP("omp for reduction(+:sumV1,sumV2)")
#endif
            for(i=0; i<communitySize ; i++)
            {
                if(eigenVectorOld[i] > 0)
                {
                    v2C[vertex[i]] = newCommunity;
                    sumV1++;
                }
                else
                    sumV2++;
            }

            /* Calculating new degree sums. */
#ifdef _OPENMP
OMP("omp for reduction(+:degreeSum1, degreeSum2) private(comm)")
#endif
            for(i=0; i<communitySize; i++)
            {
                comm = v2C[vertex[i]];
                if (comm == curCommunity)
                    degreeSum1 += G->numEdges[vertex[i]+1] 
                        - G->numEdges[vertex[i]];
                else
                    degreeSum2 += G->numEdges[vertex[i]+1] 
                        - G->numEdges[vertex[i]];
            }

            /* Calculating new modularity value. */
#ifdef _OPENMP
OMP("omp for private(u,v,degree_u, degree_v, comm) reduction(+:new_modularity)")
#endif
            for(i=0; i<communitySize; i++)
            {
                u = vertex[i];
                comm = v2C[u];
                degree_u = G->numEdges[u+1] - G->numEdges[u];
                for(j=G->numEdges[u]; j< G->numEdges[u+1]; j++)
                {
                    v = G->endV[j];
                    degree_v = G->numEdges[v+1] - G->numEdges[v];
                    if ((v2C[v] == curCommunity && comm ==  curCommunity) 
                            || (v2C[v] == newCommunity && comm == newCommunity))
                        contribution[u] -=1.0;
                    else if((v2C[v] == newCommunity && comm ==  curCommunity) 
                            || (v2C[v] == curCommunity && comm == newCommunity))
                        contribution[u] += 1.0;    
                    if(comm == v2C[v])
                        new_modularity +=1.0;
                }
                if(comm == curCommunity)
                {
                    contribution[u] += (degree_u * (degreeSum1 ))/(2.0*G->m);
                    contribution[u] 
                        -= (degree_u * (degreeSum2+degree_u ))/(2.0*G->m);
                    new_modularity  
                        += -(double)(degree_u * degreeSum1)/(double)(2.0 * G->m);
                }
                else
                {
                    contribution[u] += (degree_u * (degreeSum2 ))/(2.0*G->m);
                    contribution[u] 
                        -= (degree_u * (degreeSum1 + degree_u ))/(2.0*G->m);
                    new_modularity 
                        += -(double)(degree_u * degreeSum2)/(double)(2.0 * G->m);
                }
            }
        }

        new_modularity /= (2.0*G->m);
        if(new_modularity < modularity)
        {
            for(i=0; i<communitySize; i++)
                v2C[vertex[i]] = curCommunity;
            continue;
        }
        /* Now updating the degree Vectors */
#ifdef _OPENMP
OMP("omp parallel for private(j,comm) if (communitySize>100)")
#endif
        for(i=0; i<communitySize ; i++)
        {
            comm = v2C[vertex[i]];
            degree[vertex[i]] = 0;
            for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
            {
                if(v2C[G->endV[j]] == comm)
                    degree[vertex[i]]++;
            }
        }

        /* KL - improvement */
        if(use_improvement == 1)
        {
            max_contrib = -999999;
            maxv = -1;
            flag = 0;
            flag_counter = 0;
            while(flag == 0)
            {
                flag_counter ++;
                for(i=0; i<communitySize; i++)
                {
                    if (contribution[vertex[i]]>0 && 
                            contribution[vertex[i]] > max_contrib)
                    {
                        max_contrib = contribution[vertex[i]];
                        maxv = vertex[i];
                    }
                }
                if(maxv == -1)
                    flag = 1;
                else
                {
                    /* swap communities. */
                    if(v2C[maxv] == curCommunity)
                        v2C[maxv] = newCommunity;
                    else
                        v2C[maxv] = curCommunity;
                    /* now update the neighbours */
                    for(j=G->numEdges[maxv]; j<G->numEdges[maxv+1]; j++)
                    {
                        if(v2C[G->endV[j]] == v2C[maxv])
                            contribution[G->endV[j]] -= 1.0;    
                        /* now they are in same community */
                        else if (v2C[G->endV[j]] == newCommunity || 
                                v2C[G->endV[j]] == curCommunity)
                            contribution[G->endV[j]] += 1.0;     
                        /* not they are in different community but earlier same. */
                    }
                    degree_u = G->numEdges[maxv+1] - G->numEdges[maxv]; 
#ifdef _OPENMP
OMP("omp parallel for private(v, degree_v) if (communitySize>100)")
#endif
                    for(j=0; j<communitySize; j++)
                    {
                        v = vertex[j];
                        degree_v = G->numEdges[v+1] - G->numEdges[v];
                        if(v2C[v] == v2C[maxv])            /*same community now */
                            contribution[v] += (degree_v*degree_u)/(2.0*G->m);
                        else              /* different community earlier same. */
                            contribution[v] -= (degree_v*degree_u)/(2.0*G->m);
                    }
                }
                contribution[maxv] = max_contrib = -99999;
                maxv = -1;
            }
        }
        *numCommunities = *numCommunities + 1;
        append(Q,makeNode(curCommunity));
        append(Q,makeNode(newCommunity));
    }

}



void computeEigen(graph_t *G, double *eigenVectorOld, double *eigenVectorNew, 
        attr_id_t *v2C, attr_id_t *v2pos, attr_id_t* degree, attr_id_t
        *vertex, attr_id_t *toSplit, attr_id_t currCommunity, 
        attr_id_t communitySize, attr_id_t degreeSum) {
    attr_id_t i,j;
    attr_id_t iterCount,count, niter;
    double eigenValue,degree_u;
    double normalizedSum,ktx,mneg;
    attr_id_t numThreads;

    niter = (communitySize > 100) ? communitySize:100;
    mneg = 0.0;
    niter = 10*log(communitySize);
    count = 0;

    while(1)
    {
        iterCount = 0;
        count++;
        normalizedSum=0.0;
        srand(time(NULL));
#ifdef _OPENMP
OMP("omp parallel if(communitySize>100) ")
#endif
        {
#ifdef _OPENMP
OMP("omp for reduction(+:normalizedSum)")
#endif
            for(i=0; i<communitySize; i++)
            {
                eigenVectorOld[i] = 2.0 *((double)rand()/(double)RAND_MAX) 
                    -1.0;
                normalizedSum += eigenVectorOld[i] * eigenVectorOld[i];
            }
#ifdef _OPENMP    
OMP("omp single")
#endif
            {
                normalizedSum = sqrt(normalizedSum);
            }
#ifdef _OPENMP
OMP("omp for")
#endif
            for(i=0; i<communitySize; i++)
                eigenVectorOld[i] = eigenVectorOld[i]/normalizedSum;

        }

        while(iterCount <niter)
        {
            {
                iterCount++;
                ktx=0.0;
#ifdef _OPENMP
                numThreads = omp_get_num_threads();
#else
                numThreads = 1;
#endif
            }
#ifdef _OPENMP    
OMP("omp parallel if (communitySize>100)")
#endif
            {
#ifdef _OPENMP
OMP("omp for private(j,degree_u) reduction(+:ktx)")
#endif
                for(i=0; i<communitySize; i++)
                {
                    eigenVectorNew[i]=0.0;
                    degree_u = G->numEdges[vertex[i]+1] - G->numEdges[vertex[i]];
                    for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
                    {
                        if(v2C[G->endV[j]] == currCommunity)
                        {
                            eigenVectorNew[i] += G->dbl_weight_e[j] * eigenVectorOld[v2pos[G->endV[j]]];
                        }
                    }
                    eigenVectorNew[i] -= (((double)degree[vertex[i]]) - (double)(degree_u * degreeSum)/(double)(2.0*G->m))* eigenVectorOld[i];
                    ktx += (double) degree_u * eigenVectorOld[i];
                }
#ifdef _OPENMP
OMP("omp single")
#endif
                {
                    ktx /=(double)(2.0 * G->m);
                    normalizedSum = 0.0;
                }
#ifdef _OPENMP
OMP("omp for reduction(+:normalizedSum)")
#endif
                for(i=0; i<communitySize; i++)
                {
                    eigenVectorNew[i] -= (double)(G->numEdges[vertex[i]+1]-G->numEdges[vertex[i]]) *ktx;
                    eigenVectorNew[i] -= mneg*eigenVectorOld[i];
                    normalizedSum += eigenVectorNew[i]*eigenVectorNew[i];
                }
#ifdef _OPENMP
OMP("omp single")
#endif
                {
                    normalizedSum = sqrt(normalizedSum);
                }
#ifdef _OPENMP
OMP("omp for")
#endif
                for(i=0; i<communitySize;i++)
                {
                    eigenVectorOld[i] = eigenVectorNew[i]/normalizedSum;
                }
            }
        }
        eigenValue =0.0;
        ktx=0.0;
#ifdef _OPENMP
OMP("omp parallel if (communitySize>100)")
#endif
        {
#ifdef _OPENMP
OMP("omp for reduction(+:ktx,eigenValue) private(j,degree_u)")
#endif
            for(i=0; i<communitySize; i++)
            {
                degree_u = G->numEdges[vertex[i]+1] - G->numEdges[vertex[i]];
                for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
                {
                    if(v2C[G->endV[j]] == currCommunity)
                    {
                        eigenValue += G->dbl_weight_e[j] * eigenVectorOld[v2pos[G->endV[j]]] * eigenVectorOld[i];
                    }
                }
                eigenValue += -((double)degree[vertex[i]] - (double)(((double) degree_u *degreeSum)/(double)(2.0*G->m)))* eigenVectorOld[i] * eigenVectorOld[i];
                ktx += (double)degree_u * eigenVectorOld[i];
            }
#ifdef _OPENMP
OMP("omp single")
#endif
            {
                ktx /=(double)(2.0 * G->m);
            }
#ifdef _OPENMP
OMP("omp for reduction(-:eigenValue)")
#endif
            for(i=0; i<communitySize; i++)
                eigenValue -= ((double)(G->numEdges[vertex[i]+1]- G->numEdges[vertex[i]])) *ktx *eigenVectorOld[i];

        }    
        /* printf("The eigenValue is %f\n",eigenValue); */
        {
            if(eigenValue <0.0000001)
            {
                if(count==2)
                {
                    *toSplit = 0;
                    break;
                }
                mneg=eigenValue;
            }
            else
            {    
                *toSplit = 1;
                break;
            }
        }
    }    
}





/* Simplistic mod without klin */
void modularity_spectral_wo_klin(graph_t *G, attr_id_t *membership, attr_id_t *numCommunities)
{
    attr_id_t *v2C, *degree, *vertex, *v2pos;
    attr_id_t curCommunity=0, newCommunity,toSplit;
    attr_id_t n=G->n,sumV1,sumV2,comm,count1,count2;
    attr_id_t i,j,communitySize,degreeSum;
    list_t *Q;
    node_t *first;
    double *eigenVectorOld,*eigenVectorNew;
    attr_id_t continue_flag = 0;    

    *numCommunities = 1;

    v2C = membership;         /* v2C is a map from vertex to the 
                                 community it belongs. */
    vertex = (attr_id_t *) malloc(sizeof(attr_id_t)*n);  /* vertex is an 
                                                            array of vertices belonging to a particular community. */
    v2pos  = (attr_id_t *) malloc(sizeof(attr_id_t)*n);  
    /* v2pos is a map from vertex to its respective position in the community.*/
    degree =  (attr_id_t *) malloc(sizeof(attr_id_t)*n); 
    /* degree is a map from vertex to its degree in its respective community.*/

    eigenVectorOld=(double*)malloc(sizeof(double)*n);
    eigenVectorNew=(double*)malloc(sizeof(double)*n);
    assert(eigenVectorOld!=NULL);assert(eigenVectorNew!=NULL);

    for(i=0; i< G->n; i++)
    {
        v2C[i] = 0;
        vertex[i] = i;
        v2pos[i] = i;
        degree[i] = G->numEdges[i+1] - G->numEdges[i];
    }

    /* Making a queue. This queue will store all the 
       communities that are yet to be processed. */
    Q=(list_t*)makeList();
    append(Q, makeNode(curCommunity));

    while(Q->size > 0)
    {
        first = (node_t*) getFirst(Q);
        curCommunity = first->id;
        deleteFirst(Q);
        continue_flag = 0;

        for(i=0; i<G->n; i++)
        {
            vertex[i] = -1;
        }

        /* printf("\n\nEvaluating Community:%d\n",curCommunity); */
        degreeSum=0;
        communitySize=-1;

        /* Checking which all vertices belong to this community 
           and updating the vertex Vector accordingly. */
        /* #pragma omp parallel for shared(communitySize) 
           reduction(+:degreeSum) */
        for(i=0; i<G->n; i++)
        {
            if(v2C[i] == curCommunity)
            {
                {
                    communitySize++;
                    vertex[communitySize] = i;
                }
                v2pos[i] = communitySize;
                degreeSum += G->numEdges[i+1]-G->numEdges[i];
            }
        }
        communitySize ++;
        /* printf("community Size =%d, degree Sum =%d\n"
           ,communitySize, degreeSum); */
        if(communitySize == 1)    continue;
        computeEigen(G,eigenVectorOld, eigenVectorNew, v2C,v2pos,degree,vertex,&toSplit,curCommunity,communitySize,degreeSum);

        if(toSplit == 0)
            continue;

        newCommunity = *numCommunities;
        count1=count2=sumV1=sumV2=0;
#ifdef _OPENMP    
OMP("omp parallel for reduction(+:count1,count2)")
#endif
        for(i=0; i<communitySize; i++)
        {
            if(eigenVectorOld[i] > 0) count1++;
            else count2++;
        }
        if(count1 == 0 || count2 == 0)
        {
            continue;     /*All eigen values are of same size and 
                            hence no division is required. */
        }
#ifdef _OPENMP
OMP("omp parallel if (communitySize>100)")
#endif
        {
#ifdef _OPENMP
OMP("omp for reduction(+:sumV1,sumV2)")
#endif
            for(i=0; i<communitySize ; i++)
            {
                if(eigenVectorOld[i] > 0)
                {
                    v2C[vertex[i]] = newCommunity;
                    sumV1++;
                }
                else
                    sumV2++;
            }
            /* Now updating the degree Vectors */
#ifdef _OPENMP
OMP("omp for private(j,comm)")
#endif
            for(i=0; i<communitySize ; i++)
            {
                comm = v2C[vertex[i]];
                degree[vertex[i]] = 0;
                for(j=G->numEdges[vertex[i]]; j<G->numEdges[vertex[i]+1]; j++)
                {
                    if(v2C[G->endV[j]] == comm)
                        degree[vertex[i]]++;
                }
            }

        }            
        *numCommunities = *numCommunities + 1;
        append(Q,makeNode(curCommunity));
        append(Q,makeNode(newCommunity));


    }
}

