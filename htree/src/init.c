#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */

extern void quantile_aux(double *, double *, int *, double *, double *);
extern void quantile_R(double *, int *, double *, int *, double *, int *, double *);
extern void read_basic(double *,int *,double *, int *, int *,int *,int *,int *, int *,int *,int *,int *,int *,int *);
extern void free_basic();
extern void mtree(int *,int *,double *,double *,double *);
extern void create_history(double *, double *, int *,int *,double *,int *,int *);
extern void ctree(double *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,int *,double *,double *,int *,int *,double *,int *,int *,int *,int *,double *,int *,int *,double *,double *,double *,double *,int *,double *,int *);
extern void ctree_predict(double *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,double *,double *,int *,int *,int *,double *,int *,double *,int *,int *,double *);
extern void ctree_varimp(double *,int *,double *,int *,int *,int *,int *,int *,int *,int *,int *,double *,double *,int *,int *,int *,double *,int *,double *,int *,int *,int *,int *,double *);

static const R_CMethodDef CEntries[] = {
   {"quantile_aux",          (DL_FUNC) &quantile_aux,           5},
    {"quantile_R",            (DL_FUNC) &quantile_R,             7},
    {"read_basic",            (DL_FUNC) &read_basic,   14},
    {"free_basic",            (DL_FUNC) &free_basic,   0},
    {"mtree",                 (DL_FUNC) &mtree,   5},
    {"create_history",                 (DL_FUNC) &create_history,   7},
    {"ctree",                 (DL_FUNC) &ctree,   33},
    {"ctree_predict",         (DL_FUNC) &ctree_predict, 22},
    {"ctree_varimp",         (DL_FUNC) &ctree_varimp, 24},
    {NULL, NULL, 0}
};



void R_init_htree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
