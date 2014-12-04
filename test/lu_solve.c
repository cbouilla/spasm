#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A;
  spasm_GFp *x, *y, *b;
  int n, m, test, i, prime, result;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_triplet(stdin, 7);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  prime = A->prime;

  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));
  b = malloc(m * sizeof(spasm_GFp));

  /* test A ------------------------- with a sensible RHS ----------- */
  printf("# testing correct solution\n");
  for(i = 0; i < n; i++) {
    x[i] = rand() % prime;
  }
  for(i = 0; i < m; i++) {
    b[i] = 0;
  }
  spasm_gaxpy(A, x, b);

  result = spasm_LU_solve(A, b, x);
  if (result != SPASM_SUCCESS) {
    printf("not ok %d - LU solver [solution not found]\n", test);
    exit(1);
  }

  for(i = 0; i < m; i++) {
    y[i] = 0;
  }
  spasm_gaxpy(A, x, y);
  for(i = 0; i < m; i++) {
    if (y[i] != b[i]) {
      printf("not ok %d - LU solver [incorrect solution found]\n", test);
      exit(1);
    }
  }

  /* test B ------------------------- with a bogus RHS ----------- */
  if (n < m) {
    printf("# testing bogus solution\n");
    for(i = 0; i < m; i++) {
      b[i] = rand() % prime;
    }

    result = spasm_LU_solve(A, b, x);
    if (result == SPASM_SUCCESS) {
      printf("not ok %d - LU solver [bogus solution found]\n", test);
      exit(1);
    }

  }

  printf("ok %d - LU solver\n", test);

  spasm_csr_free(A);
  free(x);
  free(y);
  free(b);
  return 0;
}
