#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"

/** Compute the clustering coefficient for a single community.

  @param g A graph_t, either directed or undirected.  The computations
    differ for directed v. undirected graphs.
  @param membership An array mapping each vertex in g to its component.
  @param community The community of interest.

  @return Clustering coefficient.  The extreme cases without any
    triangles or open triples returns 0.
*/
double
get_single_community_clustering_coefficient (const graph_t *g,
					     const attr_id_t * membership,
					     attr_id_t community)
{
     //fprintf(stderr, "in get clustering coefficients\n");
     attr_id_t u;

     long numtri = 0;
     long numopen = 0;
     double retval;

     OMP("omp parallel for reduction(+:numopen) reduction(+:numtri)")
	  for(u=0; u<g->n; u++) {
	       const attr_id_t u_edge_start = g->numEdges[u];
	       const attr_id_t u_edge_end = g->numEdges[u+1];
	       attr_id_t u_degree = 0, j;

	       if(membership[u] != community) continue;

	       for(j=u_edge_start; j<u_edge_end; j++)
	       {
		    const attr_id_t v = g->endV[j];
		    if(membership[v] == community)
			 u_degree++;
	       }
	       for(j=u_edge_start; j<u_edge_end; j++)
	       {
		    const attr_id_t v = g->endV[j];
		    if(membership[v] == community)
		    {
			 const attr_id_t v_edge_start = g->numEdges[v];
			 const attr_id_t v_edge_end = g->numEdges[v+1];
			 attr_id_t v_degree = 0;
			 attr_id_t k = u_edge_start;
			 attr_id_t l = v_edge_start;
			 int nt = 0;

			 while(k < u_edge_end && l < v_edge_end)
			 {
			      const attr_id_t w = g->endV[k];
			      const attr_id_t x = g->endV[l];
			      if(membership[w] == community) v_degree++;
			      if(membership[w] == community && w == x) nt++;
			      if(w == x)
			      {
				   k++;
				   l++;
			      }
			      else if(w > x) l++;
			      else k++;
			 }
			 numopen += v_degree;
			 numtri += nt;
		    }
	       }
	  }

     assert (numopen >= 0);
     assert (numtri >= 0);
     if (numopen) {
	  retval = numtri/(double)numopen;
     } else if (numtri)
	  retval = HUGE_VAL;
     else
	  retval = 0.0;
     assert (retval >= 0.0);
     return retval;
}
