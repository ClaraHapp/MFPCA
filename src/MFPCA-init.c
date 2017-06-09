// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "MFPCA.h"

static const
R_CMethodDef cMethods[] = {
  {"calcCoefs",  (DL_FUNC) &calcCoefs, 4},
  {"calcImage", (DL_FUNC) &calcImage, 4},
  {"calcCoefs3D",  (DL_FUNC) &calcCoefs3D, 3},
  {"calcImage3D", (DL_FUNC) &calcImage3D, 3},
  NULL
};

void R_init_MFPCA(DllInfo* info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}