#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *U, *L;
  spasm_lu *PLUQ;
  spasm_GFp *x, *y, *u, *v, *w, *z;
  int n, m, r;
  int *p, *qinv;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  qinv = spasm_malloc(m * sizeof(int));
  p = spasm_malloc(n * sizeof(int));
  spasm_find_pivots(A, p, qinv);  

  PLUQ = spasm_PLUQ(A, p, SPASM_KEEP_L);
  
  U = PLUQ->U;
  L = PLUQ->L;
  r = U->n;

  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));
  u = malloc(n * sizeof(spasm_GFp));
  v = malloc(r * sizeof(spasm_GFp));
  w = malloc(m * sizeof(spasm_GFp));
  z = malloc(m * sizeof(spasm_GFp));

  if (spasm_is_upper_triangular(U)) {
    printf("# PLUQ = P'A ---> U is upper\n");
  } else {
    printf("not ok - PLUQ = P'A : U is not upper-triangular\n");
    exit(1);
  }

  if (spasm_is_lower_triangular(L)) {
    printf("# PLUQ = P'A ---> L is lower\n");
  } else {
    printf("not ok - PLUQ = P'A : L is not lower-triangular\n");
    exit(1);
  }

  for(int i = 0; i < n; i++) {
    spasm_vector_set(x, 0, n, 0);
    spasm_vector_set(v, 0, r, 0);
    spasm_vector_set(y, 0, m, 0);
    spasm_vector_set(w, 0, m, 0);
    spasm_vector_set(z, 0, m, 0);
    x[i] = 1;

    spasm_gaxpy(A, x, y);             // y <- x*A
    spasm_pvec(p, x, u, n);
    spasm_pvec(PLUQ->p, u, x, n);     // u <--- x.P
    spasm_gaxpy(L, x, v);             // v <--- x.(P.L)
    spasm_gaxpy(U, v, w);             // w <--- x.(P.L.U)
    spasm_pvec(PLUQ->qinv, w, z, m);  // u <--- x.(P.L.U.Q)

    for(int j = 0; j < m; j++)
      if (y[j] != z[j]) {
	      printf("not ok - PLUQ == P'A (col %d)\n", j);
	      exit(1);
      }
  }
  printf("ok - PLUQ == P'A\n");

  spasm_csr_free(A);
  spasm_free_LU(PLUQ);
  free(x);
  free(y);
  free(u);
  free(v);
  free(w);
  free(z);
  free(p);
  free(qinv);
  return 0;
}
