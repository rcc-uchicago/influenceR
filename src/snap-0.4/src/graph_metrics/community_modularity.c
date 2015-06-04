#include "graph_defs.h"
#include "graph_metrics.h"
#include "utils.h"

#if !defined(NAN)
#define NAN (0.0/0.0)
#endif

static double
get_community_modularity_dir(const graph_t *g, const attr_id_t *membership,
			     attr_id_t num_components)
{
     attr_id_t u;
     attr_id_t n, m;
     double m_inv;
     double mod;
     double *A, *E;
     attr_id_t *Lss, *Lsplus, *Lpluss;

     n = g->n;
     m = g->m;
     m_inv = 1.0/(g->m);

     Lss = calloc (3*num_components, sizeof (*Lsplus));
     assert (Lss);
     Lsplus = &Lss[num_components];
     Lpluss = &Lsplus[num_components];

     OMP("omp parallel for")
	  for(u=0; u<n; u++) {
	       const attr_id_t cid = membership[u];
	       const attr_id_t u_edge_start = g->numEdges[u];
	       const attr_id_t u_edge_end   = g->numEdges[u+1];
	       attr_id_t lss = 0, lsplus = 0;

	       if (cid >= 0) {
		    attr_id_t j;
		    for (j=u_edge_start; j<u_edge_end; j++) {
			 const attr_id_t v = g->endV[j];
			 const attr_id_t vcid = membership[v];
			 if (vcid == cid)
			      ++lss;
			 else {
			      ++lsplus;
			      if (vcid >= 0)
				   OMP("omp atomic") ++Lpluss[vcid];
			 }
		    }
		    OMP("omp atomic") Lss[cid] += lss;
		    OMP("omp atomic") Lsplus[cid] += lsplus;
	       }
	  }

     mod = 0.0;
     {
	  attr_id_t j;
	  for (j = 0; j < num_components; j++) {
	       mod += Lss[j] - (m_inv * Lsplus[j]) * Lpluss[j];
	  }
     }
     mod = mod * m_inv;

     free(Lss);

     return mod;
}

static double
get_community_modularity_undir(const graph_t *g, const attr_id_t *membership,
			       attr_id_t num_components)
{
     attr_id_t u;
     attr_id_t n, m;
     double m_inv;
     double mod;
     double *A, *E;
     attr_id_t *Lss, *Lsplus;

     n = g->n;
     m = g->m;
     m_inv = 1.0/(g->m);

     Lss = calloc (2*num_components, sizeof (*Lsplus));
     assert (Lss);
     Lsplus = &Lss[num_components];

     OMP("omp parallel for")
	  for(u=0; u<n; u++) {
	       const attr_id_t cid = membership[u];
	       const attr_id_t u_edge_start = g->numEdges[u];
	       const attr_id_t u_edge_end   = g->numEdges[u+1];
	       attr_id_t lss = 0, lsplus = 0;

	       if (cid >= 0) {
		    attr_id_t j;
		    for (j=u_edge_start; j<u_edge_end; j++) {
			 const attr_id_t v = g->endV[j];
			 const attr_id_t vcid = membership[v];
			 /* Examine only the upper triangle and diagonal. */
			 if (v < u) continue;
			 if (vcid == cid)
			      ++lss;
			 else
			      ++lsplus;
		    }
		    OMP("omp atomic") Lss[cid] += lss;
		    OMP("omp atomic") Lsplus[cid] += lsplus;
	       }
	  }

     mod = 0.0;
     {
	  attr_id_t j;
	  for (j = 0; j < num_components; j++) {
	       mod += Lss[j] - (m_inv * Lsplus[j]) * (0.25*Lsplus[j]);
	  }
     }
     mod = mod * m_inv;

     free(Lss);

     return mod;
}

/** Compute the composite modularity of the decomposition in the
  membership array.  The mechanics of these modularity computations
  are described in McCloskey and Bader, "Modularity and Graph
  Algorithms", presented at UMBC Sept. 2009.

  @param g A graph_t, either directed or undirected.  The computations
    differ for directed v. undirected graphs.
  @param membership An array mapping each vertex in g to its component.
  @param num_components The number of components, at least as large as
    the largest number in membership.

  @return The modularity, or NaN if the arguments are invalid.
*/
double
get_community_modularity(const graph_t *g, const attr_id_t *membership,
			 attr_id_t num_components)
{
     if (!g || !membership || num_components <= 0) return NAN;
     if (!g->n || !g->m) return 0.0;
     if (g->undirected)
	  return get_community_modularity_undir (g, membership, num_components);
     else
	  return get_community_modularity_dir (g, membership, num_components);
}

static double
get_single_community_modularity_undir(const graph_t *g, const attr_id_t *membership,
				      attr_id_t which_component)
{
     attr_id_t n, m;
     double mod;
     attr_id_t Ls = 0, Vols = 0, Xs, u;

     n = g->n;
     m = g->m/2;

     OMP("omp parallel for reduction(+:Ls) reduction(+:Vols)")
	  for(u=0; u<n; u++) {
	       const attr_id_t u_edge_start = g->numEdges[u];
	       const attr_id_t u_edge_end   = g->numEdges[u+1];
	       const attr_id_t u_degree     = u_edge_end - u_edge_start;
	       attr_id_t j;

	       if (membership[u] == which_component) {
		    for(j=u_edge_start; j<u_edge_end; j++) {
			 const attr_id_t v = g->endV[j];
			 if (membership[v] == which_component)
			      ++Ls; /* double-counted */
		    }
		    Vols += u_degree;
	       }
	  }

     assert (!(Ls % 2));
     Ls /= 2;
     Xs = Vols - Ls;
     mod = (Ls - (Xs*(double)Xs)/(4.0*m))/m;

     return mod;
}

static double
get_single_community_modularity_dir(const graph_t *g, const attr_id_t *membership,
				    attr_id_t which_component)
{
     attr_id_t n, m;
     double mod;
     attr_id_t Ls = 0, Vols = 0, Xs, u;

     n = g->n;
     m = g->m;

     OMP("omp parallel for reduction(+:Ls) reduction(+:Vols)")
	  for(u=0; u<n; u++) {
	       const attr_id_t u_edge_start = g->numEdges[u];
	       const attr_id_t u_edge_end   = g->numEdges[u+1];
	       const attr_id_t u_degree     = u_edge_end - u_edge_start;
	       attr_id_t j;

	       if (membership[u] == which_component) {
		    for(j=u_edge_start; j<u_edge_end; j++) {
			 const attr_id_t v = g->endV[j];
			 if (membership[v] == which_component)
			      ++Ls; /* double-counted */
		    }
		    Vols += u_degree;
	       }
	  }

     assert (!(Ls % 2));
     Ls /= 2;
     Xs = Vols - Ls;
     mod = (Ls - (Xs*(double)Xs)/m)/m;

     return mod;
}

/** Compute the composite modularity of a single component according
  to the decomposition in the membership array.  The
  get_community_modularity() routine is more efficient than computing
  each modularity contribution separately. The mechanics of these
  modularity computations are described in McCloskey and Bader,
  "Modularity and Graph Algorithms", presented at UMBC Sept. 2009.

  @param g A graph_t, either directed or undirected.  The computations
    differ for directed v. undirected graphs.
  @param membership An array mapping each vertex in g to its component.
  @param which_component The component of interest.

  @return The modularity, or NaN if the arguments are invalid.
*/
double
get_single_community_modularity(const graph_t *g, const attr_id_t *membership,
				attr_id_t which_component)
{
     if (!g || !membership) return NAN;
     if (g->undirected)
	  return get_single_community_modularity_undir (g, membership, which_component);
     else
	  return get_single_community_modularity_dir (g, membership, which_component);
}
