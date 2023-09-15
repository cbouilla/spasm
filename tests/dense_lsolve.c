#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *G;
  int i, n, m, prime;
  spasm_ZZp *x, *b, *y;
  int test = 0;
  
  T = spasm_load_sms(stdin, 32003);
  G = spasm_compress(T);
  spasm_triplet_free(T);

  n = G->n;
  m = G->m;
  prime = G->field.p;
  assert(n >= m);
  assert(spasm_is_lower_triangular(G));

  x = malloc(n * sizeof(spasm_ZZp));
  b = malloc(m * sizeof(spasm_ZZp));
  y = malloc(m * sizeof(spasm_ZZp));

  for(i = 0; i < m; i++) {
    b[i] = spasm_ZZp_init(&G->field, rand());
    y[i] = b[i];
  }

  spasm_dense_back_solve(G, y, x, SPASM_IDENTITY_PERMUTATION);

  for(i = 0; i < m; i++) {
    y[i] = 0;
  }
  spasm_xApy(x, G, y);
  for(i = 0; i < m; i++) {
    if (y[i] != b[i]) {
      printf("not ok %d - dense back-substitution triangular solver [incorrect solution found]\n", test);
      exit(1);
    }
  }

  printf("ok %d -  dense back-substitution triangular solver\n", test);

  spasm_csr_free(G);
  free(x);
  free(y);
  free(b);
  return 0;
}
