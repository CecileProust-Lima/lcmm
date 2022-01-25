#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "lcmm.h"

static R_FortranMethodDef FortRout[] = {
  {"hetmixlin", (DL_FUNC) &F77_SUB(hetmixlin), 35},
  {"hetmixcont", (DL_FUNC) &F77_SUB(hetmixcont), 47},
  {"hetmixord", (DL_FUNC) &F77_SUB(hetmixord), 39},
  {"jointhet", (DL_FUNC) &F77_SUB(jointhet), 59},
  {"predictcont", (DL_FUNC) &F77_SUB(predictcont), 21},
  {"predictmult", (DL_FUNC) &F77_SUB(predictmult), 25},
  {"predictcondmult", (DL_FUNC) &F77_SUB(predictcondmult), 17},
  {"cvpl", (DL_FUNC) &F77_SUB(cvpl), 41},
  {"postprob2", (DL_FUNC) &F77_SUB(postprob2), 40},
  {"calculustransfo", (DL_FUNC) &F77_SUB(calculustransfo), 16},
  {"mpjhet", (DL_FUNC) &F77_SUB(mpjhet), 64},
  {"hetmixmult", (DL_FUNC) &F77_SUB(hetmixmult), 59},
  {"iteminfo", (DL_FUNC) &F77_SUB(iteminfo), 15},
  {NULL, NULL, 0}
};


void R_init_lcmm(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
