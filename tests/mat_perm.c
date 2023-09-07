#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "spasm.h"

int main(int argc, char **argv) {
  int n, m, i, j;
  spasm_triplet *T;
  spasm *A, *B;
  int *p, *q;
  spasm_GFp *x, *y, *u, *v, *w;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  // generate random row & col permutation
  p = spasm_random_permutation(n);
  q = spasm_random_permutation(m);

  B = spasm_permute(A, p, q, true);

  x = spasm_malloc(n * sizeof(spasm_GFp));
  y = spasm_malloc(m * sizeof(spasm_GFp));
  u = spasm_malloc(n * sizeof(spasm_GFp));
  v = spasm_malloc(m * sizeof(spasm_GFp));
  w = spasm_malloc(m * sizeof(spasm_GFp));

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      x[j] = 0;
    }
    for(j = 0; j < m; j++) {
      y[j] = 0;
      v[j] = 0;
    }
    x[i] = 1;

    spasm_xApy(x, A, y); // y <--- x.A

    spasm_pvec(p, x, u, n); // u <--- x.P
    spasm_xApy(u, B, v); // v <--- (x.P).B
    spasm_pvec(q, v, w, m); // w <--- ((x.P).B).Q

    for(j = 0; j < m; j++) {
      if (y[j] != w[j]) {
	printf("not ok - P*A*Q \n");
	exit(1);
      }
    }
  }
  printf("ok - P*A*Q \n");

  free(p);
  free(q);
  free(x);
  free(y);
  free(u);
  free(v);
  free(w);
  spasm_csr_free(A);
  spasm_csr_free(B);
}
