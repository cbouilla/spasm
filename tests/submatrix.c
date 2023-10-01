#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *S;
  int n, m;

  T = spasm_load_sms(stdin, 46337, NULL);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  assert(n > 10);
  assert(m > 10);

  S = spasm_submatrix(A, 5, 10, 5, 10, true);

  spasm_save_csr(stdout, S);
  spasm_csr_free(S);
  spasm_csr_free(A);
}
