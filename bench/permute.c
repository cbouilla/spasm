#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *PAP;
  FILE *f;
  int n;
  int *p, *pinv;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  assert(n == A->m);
  assert(argc > 1);
  
  f = fopen(argv[1], "r");
  p = spasm_load_permutation(f, n);
  fclose(f);
  pinv = spasm_pinv(p, n);
  PAP = spasm_permute(A, pinv, p, SPASM_WITH_NUMERICAL_VALUES);
  
  spasm_save_csr(stdout, PAP);

  spasm_csr_free(PAP);
  spasm_csr_free(A);
  free(p);
  free(pinv);
  return 0;
}
