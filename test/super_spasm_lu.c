#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *M;
  super_spasm *A;
  spasm_lu *LU;
  spasm_GFp *x, *y, *u, *v;
  int n, m, test, i, j;

  assert(argc >1);
  test = atoi(argv[1]);

  T = spasm_load_sms(stdin, 42013);
  A = spasm_malloc(sizeof(super_spasm));
  A->M = M = spasm_compress(T);
  A->n = M->n + 2;
  A->p = spasm_malloc(M->n * sizeof(int));
  n = M->n;
  m = M->m;
  spasm_triplet_free(T);

  for(i = 0; i < n; i++){
    A->p[i] = i+1;
  }

  //Get workspace
  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));
  u = malloc(n * sizeof(spasm_GFp));
  v = malloc(m * sizeof(spasm_GFp));

  LU = super_spasm_LU(A, 0, SPASM_KEEP_L);

  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      x[j] = 0;
      u[j] = 0;
    }
    for(j = 0; j < m; j++) {
      y[j] = 0;
      v[j] = 0;
    }
    x[i] = 1;

    spasm_gaxpy(M, x, y);     // y <- x*M
    spasm_gaxpy(LU->L, x, u); // u <- x*L
    spasm_gaxpy(LU->U, u, v); // v <- (x*L)*U

    for(j = 0; j < m; j++) {
      if(y[j] != v[j]) {
	printf("not ok %d - L*U == M (col %d)\n", test, j);
	exit(0);
      }
    }
  }

  printf("ok %d - L*U == A\n", test);

  super_spasm_free(A);
  spasm_free_LU(LU);
  free(x);
  free(y);
  free(u);
  free(v);
  return 0;

}
