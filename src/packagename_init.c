#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP process_sparse_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP snap_betweenness_R(SEXP, SEXP, SEXP);
extern SEXP snap_bridging_R(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP snap_keyplayer_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"process_sparse_R",   (DL_FUNC) &process_sparse_R,   6},
    {"snap_betweenness_R", (DL_FUNC) &snap_betweenness_R, 3},
    {"snap_bridging_R",    (DL_FUNC) &snap_bridging_R,    5},
    {"snap_keyplayer_R",   (DL_FUNC) &snap_keyplayer_R,   9},
    {NULL, NULL, 0}
};

void R_init_influenceR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
