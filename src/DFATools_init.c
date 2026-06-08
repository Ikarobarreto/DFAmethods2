#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void rdcca(void *, void *, void *, void *, void *);
extern void rdfa(void *, void *, void *, void *);
extern void rdcca_box(void *, void *, void *, void *, void *, void *);
extern void rdfa_box(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"rdcca",     (DL_FUNC) &rdcca,     5},
  {"rdfa",      (DL_FUNC) &rdfa,      4},
  {"rdcca_box", (DL_FUNC) &rdcca_box, 6},
  {"rdfa_box",  (DL_FUNC) &rdfa_box,  5},
  {NULL, NULL, 0}
};

void R_init_DFATools(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
