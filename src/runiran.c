#include <R.h>

void F77_SUB(getrand)(void) {
  GetRNGstate();
}

void F77_SUB(putrand)(void) {
  PutRNGstate();
}

double F77_SUB(runiran)(void) {
  return unif_rand();
}
