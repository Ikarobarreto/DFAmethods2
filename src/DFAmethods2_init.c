#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void rdcca(void *, void *, void *, void *, void *);
extern void rdfa(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"rdcca", (DL_FUNC) &rdcca, 5},
  {"rdfa",  (DL_FUNC) &rdfa,  4},
  {NULL, NULL, 0}
};

void R_init_DFAmethods2(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
