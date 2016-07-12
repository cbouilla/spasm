/* indent -nfbs -i2 -nip -npsl -di0 -nut rank_gplu.c */
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include "spasm.h"

/** this program demonstrate a naive dense rank computation */

int main() {
  int r, i, j, px, n, m, *Aj, *Ap, prime, *q;
  spasm_triplet *T;
  spasm *A;
  spasm_dense_lu *LU;
  spasm_GFp *x, *y, *Ax;

  prime = 42013;
  
  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;

  LU = spasm_dense_LU_alloc(m, prime);
  q = spasm_malloc(m * sizeof(int));
  x = spasm_malloc(m * sizeof(spasm_GFp));
  y = spasm_malloc(m * sizeof(spasm_GFp));
  for(j = 0; j < m; j++) {
    q[j] = j;
  }

  r = 0;
  for(i = 0; i < n; i++) {
    spasm_vector_zero(x, m);
    for(px = Ap[i]; px < Ap[i+1]; px++) {
      x[ Aj[px] ] = Ax[px];
    }
    for(j = 0; j < m; j++) {
      y[j] = x[ q[j] ];
    }
    r += spasm_dense_LU_process(LU, y, q);

    fprintf(stderr, "\rrow %d/%d, rank >= %d", i+1, n, r);
    fflush(stderr);
  }

  fprintf(stderr, "\nFinal rank = %d\n", r);

  free(q);
  free(x);
  spasm_dense_LU_free(LU);
  spasm_csr_free(A);
  return 0;
}
