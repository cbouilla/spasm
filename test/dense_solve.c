#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *G;
  int i, n, m, test, prime, result;
  spasm_GFp *x, *b, *y;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_triplet(stdin, 65537);
  G = spasm_compress(T);
  spasm_triplet_free(T);

  n = G->n;
  m = G->m;
  prime = G->prime;
  assert(n <= m);

  x = malloc(m * sizeof(spasm_GFp));
  b = malloc(m * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));

  /* test A ------------------------- with a sensible RHS ----------- */
  printf("# testing correct solution\n");

  for(i = 0; i < m; i++) {
    x[i] = rand() % prime;
    b[i] = 0;
    y[i] = 0;
  }
  spasm_gaxpy(G, x, b);
  for(i = 0; i < m; i++) {
    x[i] = b[i];
  }
  // x=b is a random RHS

  result = spasm_dense_forward_solve(G, x, SPASM_IDENTITY_PERMUTATION);
  if (result != SPASM_SUCCESS) {
    printf("not ok %d - dense forward-substitution triangular solver [solution not found]\n", test);
    exit(1);
  }

  spasm_gaxpy(G, x, y);
  for(i = 0; i < n; i++) {
    if (y[i] != b[i]) {
      printf("not ok %d - dense forward-substitution triangular solver [incorrect solution found]\n", test);
      exit(1);
    }
  }

  /* test B ------------------------- with a bogus RHS ----------- */
  if (n < m) {
    printf("# testing bogus solution\n");
    for(i = 0; i < m; i++) {
      x[i] = rand() % prime;
    }

    result = spasm_dense_forward_solve(G, x, SPASM_IDENTITY_PERMUTATION);
    if (result == SPASM_SUCCESS) {
      printf("not ok %d - dense forward-substitution triangular solver [bogus solution found]\n", test);
      exit(1);
    }
  }

  printf("ok %d -  dense forward-substitution triangular solver\n", test);

  spasm_csr_free(G);
  free(x);
  free(y);
  free(b);
  return 0;
}
