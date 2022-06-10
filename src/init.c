#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "lcmm.h"

static R_FortranMethodDef FortRout[] = {
  {"predictcont", (DL_FUNC) &F77_SUB(predictcont), 21},
  {"predictmult", (DL_FUNC) &F77_SUB(predictmult), 25},
  {"predictcondmult", (DL_FUNC) &F77_SUB(predictcondmult), 17},
  {"cvpl", (DL_FUNC) &F77_SUB(cvpl), 41},
  {"postprob2", (DL_FUNC) &F77_SUB(postprob2), 40},
  {"calculustransfo", (DL_FUNC) &F77_SUB(calculustransfo), 16},
  {"iteminfo", (DL_FUNC) &F77_SUB(iteminfo), 15},
  {"loglikmultlcmm", (DL_FUNC) &F77_SUB(loglikmultlcmm), 53},
  {"loglikjointlcmm", (DL_FUNC) &F77_SUB(loglikjointlcmm), 60},
  {"loglikmpjlcmm", (DL_FUNC) &F77_SUB(loglikmpjlcmm), 62},
  {"logliklcmmcont", (DL_FUNC) &F77_SUB(logliklcmmcont), 42},
  {"logliklcmmord", (DL_FUNC) &F77_SUB(logliklcmmord), 35},
  {"loglikhlme", (DL_FUNC) &F77_SUB(loglikhlme), 30},
  {NULL, NULL, 0}
};


void R_init_lcmm(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
