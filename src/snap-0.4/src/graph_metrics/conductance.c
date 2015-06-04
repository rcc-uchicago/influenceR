#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"

/** Compute the conductance of a single community.  Note that the
  definition is identical for directed and undirected graphs.

  @param g A graph_t, either directed or undirected.  The computations
    differ for directed v. undirected graphs.
  @param membership An array mapping each vertex in g to its component.
  @param component The component of interest.

  @return Conductance, a volume-weighted cutsize measure.  The extreme
    cases, where the community induces no cut, return HUGE_VAL.
*/
double 
get_single_community_conductance(const graph_t *g, const attr_id_t * membership,
				 attr_id_t component)
{
     attr_id_t u;
     attr_id_t degree_in = 0, degree_out = 0;
     attr_id_t cross_edges = 0;

     OMP("omp parallel for reduction(+:degree_in) reduction(+:degree_out) \
	 reduction(+:cross_edges)")
	  for(u=0; u<g->n; u++) {
	       const attr_id_t u_edge_start = g->numEdges[u];
	       const attr_id_t u_edge_end = g->numEdges[u+1];
	       const attr_id_t u_degree = u_edge_end - u_edge_start;

	       if(membership[u] == component)
	       {
		    attr_id_t j;
		    degree_in += u_degree;	

		    for(j=u_edge_start; j<u_edge_end; j++)
		    {
			 attr_id_t v = g->endV[j];
			 if(membership[v] != component)
			      cross_edges++;
		    }
	       }
	       else
	       {
		    degree_out += u_degree;	
	       }
	  }

     if(degree_in < degree_out)
     {
	  if (cross_edges)
	       return cross_edges / (double) degree_in;
	  else return HUGE_VAL;
     }
     else 
     {
	  if (cross_edges)
	       return cross_edges / (double) degree_out;
	  else return HUGE_VAL;
     }
}
