#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *L, *B;
  int i, j, n, m, test, top, *xi;
  spasm_GFp *x, *y;
  test = 0;

  // load matrix
  T = spasm_load_sms(stdin, 32003);
  L = spasm_compress(T);
  spasm_triplet_free(T);
  n = L->n;
  m = L->m;

  assert( m < n ); // lower-trapezoidal

  // load RHS
  T = spasm_triplet_alloc(1, n, n, 32003, true);
  for(j = m; j < n; j++) {
    spasm_add_entry(T, 0, j, j);
  }
  B = spasm_compress(T);
  spasm_triplet_free(T);

  xi = malloc(3*n * sizeof(int));
  spasm_vector_zero(xi, 3*n);

  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(n * sizeof(spasm_GFp));
  spasm_vector_zero(x, n);
  spasm_vector_zero(y, n);

  top = spasm_sparse_backward_solve(L, B, 0, xi, x, SPASM_IDENTITY_PERMUTATION, 0);

  spasm_gaxpy(L, x, y);
  for(j = m; j < n; j++) {
    y[j] = x[j];
  }
  spasm_scatter(B->j, B->x, B->p[0], B->p[1], B->prime - 1, y, B->prime);

  for(i = 0; i < n; i++) {
    if (y[i] != 0) {
      printf("not ok - sparse triangular L-solve (index %d, n=%d, m=%d\n", i, n, m);
      exit(1);
    }
  }

  printf("ok - sparse triangular L-solve\n");

  spasm_csr_free(L);
  spasm_csr_free(B);
  free(xi);
  free(x);
  free(y);
  exit(EXIT_SUCCESS);
}
