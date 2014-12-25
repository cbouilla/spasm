#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *L, *B;
  int i, n, m, test, top, *xi;
  spasm_GFp *x, *y;

  assert(argc == 2);
  test = atoi(argv[1]);

  // load matrix
  T = spasm_load_sms(stdin, 32003);
  L = spasm_compress(T);
  spasm_triplet_free(T);
  n = L->n;
  m = L->m;

  assert( n >= m ); // lower-trapezoidal

  // load RHS
  T = spasm_triplet_alloc(1, m, 10, 32003, true);
  spasm_add_entry(T, 0, 0, 1);
  spasm_add_entry(T, 0, m / 2, 2);
  spasm_add_entry(T, 0, m - 1, 3);
  B = spasm_compress(T);
  spasm_triplet_free(T);

  xi = malloc(3*m * sizeof(int));
  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));
  spasm_vector_zero(x, n);
  spasm_vector_zero(y, m);

  top = spasm_sparse_backward_solve(L, B, 0, xi, x, SPASM_IDENTITY_PERMUTATION);

  spasm_gaxpy(L, x, y);
  spasm_scatter(B->j, B->x, B->p[0], B->p[1], B->prime - 1, y, B->prime);

  for(i = 0; i < m; i++) {
    if (y[i] != 0) {
      printf("not ok %d - sparse triangular L-solve\n", test);
      exit(0);
    }
  }

  printf("ok %d - sparse triangular L-solve\n", test);

  spasm_csr_free(L);
  spasm_csr_free(B);
  free(xi);
  free(x);
  free(y);
  return 0;
}
