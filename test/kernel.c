#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *K;
  spasm_GFp *x, *y, *Kx;
  int i, j, n, m, k, test, p;
  int *Kp, *Kj;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  K = spasm_kernel(A);
  k = K->n;
  Kp = K->p;
  Kj = K->j;
  Kx = K->x;

  assert( K->m == A->n );

  x = spasm_malloc(n * sizeof(spasm_GFp));
  y = spasm_malloc(m * sizeof(spasm_GFp));

  /* test that they are really kernel vectors */
  for(i = 0; i < k; i++) {
    printf("# testing vector %d\n", i);
    spasm_vector_zero(x, n);
    spasm_vector_zero(y, m);

    /* scatter K[i] into x */
    for(p = Kp[i]; p < Kp[i + 1]; p++) {
      x[ Kj[p] ] = Kx[p];
    }

    /* y <-- x.A */
    spasm_gaxpy(A, x, y);

    for(j = 0; j < m; j++) {
      if (y[j] != 0) {
      printf("not ok %d - vector not in kernel\n", test);
      exit(0);
      }
    }
  }

  printf("ok %d - (left-)kernel basis\n", test);

  spasm_csr_free(A);
  spasm_csr_free(K);
  free(x);
  free(y);
  return 0;
}
