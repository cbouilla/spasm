#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A;
  spasm_lu *LU;
  spasm_ZZp *w, *x, *y, *u, *v;
  int n, m, i, j, *p, *qinv;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  w = spasm_malloc(n * sizeof(spasm_ZZp));
  x = spasm_malloc(n * sizeof(spasm_ZZp));
  y = spasm_malloc(m * sizeof(spasm_ZZp));
  u = spasm_malloc(n * sizeof(spasm_ZZp));
  v = spasm_malloc(m * sizeof(spasm_ZZp));
  
  p = spasm_malloc(n * sizeof(int));
  qinv = spasm_malloc(m * sizeof(int));
  spasm_find_pivots(A, p, qinv);
  LU = spasm_LU(A, p, SPASM_KEEP_L);

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      x[j] = 0;
      u[j] = 0;
    }
    for(j = 0; j < m; j++) {
      y[j] = 0;
      v[j] = 0;
    }
    x[i] = 1;

    spasm_ipvec(p, x, w, n);
    spasm_gaxpy(A, w, y);     // y <- x*P*A
    spasm_gaxpy(LU->L, x, u); // u <- x*L
    spasm_gaxpy(LU->U, u, v); // v <- (x*L)*U

    for(j = 0; j < m; j++) {
      if (y[j] != v[j]) {
	      printf("not ok - L*U == P*A (col %d)\n", j);
	      exit(1);
      }
    }
  }

  printf("ok - L*U == P*A\n");

  spasm_csr_free(A);
  spasm_free_LU(LU);
  free(x);
  free(y);
  free(u);
  free(v);
  free(p);
  free(qinv);
  return 0;
}
