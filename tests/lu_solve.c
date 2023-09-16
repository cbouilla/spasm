#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A;
  spasm_ZZp *x, *y, *b;
  int n, m, i, prime, result;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  prime = A->field.p;

  x = spasm_malloc(n * sizeof(spasm_ZZp));
  y = spasm_malloc(m * sizeof(spasm_ZZp));
  b = spasm_malloc(m * sizeof(spasm_ZZp));

  /* test A ------------------------- with a sensible RHS ----------- */
  printf("# testing correct solution\n");
  for (int i = 0; i < n; i++)
    x[i] = spasm_ZZp_init(A->field, rand());
  
  for (int i = 0; i < m; i++)
    b[i] = 0;
  spasm_gaxpy(A, x, b);

  for (int i = 0; i < n; i++)
    x[i] = 0;

  result = spasm_LU_solve(A, b, x);
  if (result != SPASM_SUCCESS) {
    printf("not ok - LU solver [solution not found]\n");
    exit(1);
  }

  for (int i = 0; i < m; i++)
    y[i] = 0;
  
  spasm_gaxpy(A, x, y);
  for (int i = 0; i < m; i++)
    if (y[i] != b[i]) {
      printf("not ok - LU solver [incorrect solution found]\n");
      exit(1);
    }

  /* test B ------------------------- with a bogus RHS ----------- */
  if (n < m) {
    printf("# testing bogus solution\n");
    for(i = 0; i < m; i++)
      b[i] = spasm_ZZp_init(A->field, rand());
  

    result = spasm_LU_solve(A, b, x);
    if (result == SPASM_SUCCESS) {
      printf("not ok - LU solver [bogus solution found]\n");
      exit(1);
    }

  }

  printf("ok - LU solver\n");

  spasm_csr_free(A);
  free(x);
  free(y);
  free(b);
  return 0;
}
