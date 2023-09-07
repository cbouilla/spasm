#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *G;
  int i, j, n, m, prime, result;
  spasm_GFp *x, *b, *y;

  T = spasm_load_sms(stdin, 32003);
  G = spasm_compress(T);
  spasm_triplet_free(T);

  n = G->n;
  m = G->m;
  prime = G->prime;
  assert(n <= m);

  x = malloc(n * sizeof(spasm_GFp));
  b = malloc(m * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));

  /* test A ------------------------- with a sensible RHS ----------- */
  printf("# testing correct solution\n");

  for(i = 0; i < n; i++) {
    x[i] = rand() % prime;
  }
  for(j = 0; j < m; j++) {
    b[j] = 0;
  }
  spasm_xApy(x, G, b);
  for(j = 0; j < m; j++) {
    y[j] = b[j];
  }

  result = spasm_dense_forward_solve(G, y, x, SPASM_IDENTITY_PERMUTATION);
  if (result != SPASM_SUCCESS) {
    printf("not ok - dense forward-substitution triangular solver [solution not found]\n");
    exit(1);
  }

  for(j = 0; j < m; j++) {
    y[j] = 0;
  }
  spasm_xApy(x, G, y);
  for(i = 0; i < m; i++) {
    if (y[i] != b[i]) {
      printf("not ok - dense forward-substitution triangular solver [incorrect solution found]\n");
      exit(1);
    }
  }

  /* test B ------------------------- with a bogus RHS ----------- */
  if (n < m) {
    printf("# testing bogus solution\n");
    for(j = 0; j < m; j++) {
      b[j] = rand() % prime;
    }

    result = spasm_dense_forward_solve(G, b, x, SPASM_IDENTITY_PERMUTATION);
    if (result == SPASM_SUCCESS) {
      printf("not ok - dense forward-substitution triangular solver [bogus solution found]\n");
      exit(1);
    }
  }

  printf("ok -  dense forward-substitution triangular solver\n");

  spasm_csr_free(G);
  free(x);
  free(y);
  free(b);
  exit(EXIT_SUCCESS);
}
