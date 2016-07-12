/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_lu.c  */
#include <assert.h>
#include "spasm.h"

spasm_dense_lu *spasm_dense_LU_alloc(int m, int prime) {
  spasm_dense_lu *R;

  R = spasm_malloc(sizeof(spasm_dense_lu));
  R->m = m;
  R->n = 0;
  R->prime = prime;
  R->nmax = 10;
  R->x = spasm_malloc(sizeof(spasm_GFp) * R->nmax * R->m);
  return R;
}

void spasm_dense_LU_free(spasm_dense_lu * A) {
  free(A->x);
  free(A);
}

void spasm_dense_LU_realloc(spasm_dense_lu *A, int n) {
  A->x = spasm_realloc(A->x, n * A->m * sizeof(spasm_GFp));
  A->nmax = n;
}

void spasm_dense_LU_grow(spasm_dense_lu *A, const spasm_GFp *y) {
  int n, m, j, nmax;
  spasm_GFp *x;

  n = A->n;
  m = A->m;
  nmax = A->nmax;
  if (n == nmax) {
    spasm_dense_LU_realloc(A, 2*n);
  }
  x = A->x;
  for(j = 0; j < m; j++) {
    x[n*m + j] = y[j];
  }
  A->n = n+1;
}

void spasm_dense_LU_eliminate(const spasm_dense_lu *A, spasm_GFp *y) {
  int i, j, n, m, p;
  spasm_GFp *x;

  x = A->x;
  n = A->n;
  m = A->m;
  p = A->prime;

  for (i = 0; i < n; i++) {
    for(j = i + 1; j < m; j++) {
      y[j] += (p - y[i]) * x[i*m + j];
      y[j] = y[j] % p;
    }
    y[i] = 0;
  }
}

void spasm_dense_LU_swap_column(spasm_dense_lu *A, int a, int b) {
  int i, n, m;
  spasm_GFp *x, tmp;

  x = A->x;
  n = A->n;
  m = A->m;
  
  for (i = 0; i < n; i++) {
    tmp = x[i*m + a];
    x[i*m + a] = x[i*m + b];
    x[i*m + b] = tmp;
  }
}


/* if y belongs to the linear span of U, return 0. Else update U and return 1. */
int spasm_dense_LU_process(spasm_dense_lu *A, spasm_GFp *y, int *q) {
  int j, k, n, m, p;
  spasm_GFp beta, *x;


  m = A->m;
  n = A->n;
  p = A->prime;
  x = A->x;

  spasm_dense_LU_eliminate(A, y);

  /* locate potential pivot */
  for(k = n; k < m; k++) {
    if (y[k]) {
      break;
    }
  }

  if (k == m) {
    return 0;
  }

  beta = spasm_GFp_inverse(y[k], p);
  for(j = n; j < m; j++) {
    y[j] *= beta;
    y[j] = y[j] % p;
  }

  spasm_dense_LU_grow(A, y);
  
  /* permute columns if necessary */
  if (k != n) {
    //assert(k > n);
    spasm_dense_LU_swap_column(A, k, n);
    spasm_swap(q, k, n);
  }

  return 1;
}
