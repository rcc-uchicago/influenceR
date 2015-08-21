/*
 process_sparce.c: Speed up computations on sparse matrix for constraint calculation.
 
 AUTHOR: Simon Jacobs <sdjacobs@uchicago.edu>
 LICENSE: GPLv2
*/
 

#include <R.h>
#include <Rinternals.h>

double * process_sparse(int *I, int *J, double *X, double *Ai,  double *deg, int n, double *out)
{

  for(int p = 0; p < n; p++) {
    int j = J[p];
    int i = I[p];
    out[p] = X[p] * Ai[j] * Ai[i] * deg[j];
    //out[p] = (out[p] == 0 ? 0 : 1/out[p]);
  }

  return out;
}
    
SEXP process_sparse_R(SEXP sI, SEXP sJ, SEXP sX, SEXP sAi, SEXP sdeg, SEXP sn)
{
  
  int n = INTEGER(sn)[0];
  
  SEXP sres = PROTECT(allocVector(REALSXP, n));

  /*Coerce inputs*/
  int *I = INTEGER(sI),
    *J = INTEGER(sJ);
  
  double *X = REAL(sX),
    *Ai = REAL(sAi),
    *deg = REAL(sdeg),
    *res = REAL(sres);
  
  process_sparse(I, J, X, Ai, deg, n, res);
  
  UNPROTECT(1);
 
  return sres;
}




