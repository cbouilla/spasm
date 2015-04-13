#include <assert.h>
#include "spasm.h"


/* 
 * Find x such as x * L = ek.
 *
 * The solution is scattered in x, its pattern is given by xi.
 * The return value nz is the number of non-zero entries of x. 
 */
int spasm_linverse(const spasm *L, int k, spasm_GFp *x, int *xi, const int *pinv ) {
  spasm *I;
  int n, prime, i, top, nz=0;

  n = L->n; // number of rows of L
  prime = L->prime; 

  I = spasm_csr_alloc(n, n, n, prime, 1);

  for(i=0; i<n; i++) {
    I->p[i]=i;
    I->j[i]=i;
    I->x[i]=1;
  }
  I->p[n]=n; // Identity Matrix

  top = spasm_sparse_backward_solve(L, I, k, xi, x, pinv);

  for (i = top; i < n; i++) {
    if (x[i] != 0) nz++;
  }

  spasm_csr_free(I);

  return nz;
}


/*
 * Calculate y = x * M, where x and M are sparse.
 *
 * The result is scattered in y, its pattern is given by yi. 
 * The return value nz is the number of non-zero entries in y.
 */
int spasm_sparse_vector_matrix_multiply(const spasm *M, const spasm_GFp *x, const int *xi, int xnz, spasm_GFp *y, int *yi) {
  int p, j, i=0, m, nz=0, Mnz=spasm_nnz(M), prime, *Mp, *Mj, *w, *wi;
  spasm_GFp *Mx;

  // check inputs
  if (x == NULL || Mnz == 0) {
    y = NULL;
    yi = NULL;
    return 0;
  }

  m = M->m;
  Mp = M->p;
  Mj = M->j;
  Mx = M->x;
  prime = M->prime;

  // get workspace, initializing w
  w = spasm_malloc(m * sizeof(int));
  for (j = 0; j < m; j++) {
    w[j] = -1;
  }
  wi = spasm_calloc(m, sizeof(int));

  // scatter y and find its pattern
  for (j = 0; j < xnz; j++) {
    p = xi[j];
    nz = spasm_scatter_and_pattern(Mj, Mx, Mp[p], Mp[p+1], x[p], y, wi, w, nz, prime); 
  }

  // reallocate yi
  yi = spasm_realloc(yi, nz*sizeof(int));

  for (j = 0; j < xnz; j++) {
    if (w[j] >= 0) {
      yi[i] = j;
      i++;
    }
  }

  // free workspace
  free(w);
  free(wi);
  free(Mp);
  free(Mj);
  free(Mx);
  return nz;
}


/*
 * Compress a scattered vector, given its pattern.
 */
void spasm_compress_vector(const spasm_GFp *scat, const int *xi, const int nz, spasm_GFp *comp) {
  int i;

  comp = spasm_realloc(comp, nz*sizeof(int));
  for (i = 0; i < nz; i++) {
    comp[i] = scat[xi[i]];
  }
}


/*
 * Given a matrix M and a lower_triangular matrix L, find y such as
 * y = x * M where x is solution of x * L = ek.
 *
 * the return value nz is the number of non-zero entries in y.
 */

int spasm_inverse_and_sparse_vect_multiplication(const spasm *L, const spasm *M, int k, spasm_GFp *y, int *yi, const int *pinv) {
  spasm_GFp *x, *yscat;
  int *xi, xnz, nz, Lm, m;

  Lm = L->m;
  m = M->m;

  //get workspace
  x = spasm_calloc(Lm, sizeof(spasm_GFp));
  xi = spasm_calloc(Lm, sizeof(int));
  yscat = spasm_calloc(m, sizeof(spasm_GFp));

  //solve x * L = ek
  xnz = spasm_linverse(L, k, x, xi, pinv);

  //calculate yscat = x * M
  nz = spasm_sparse_vector_matrix_multiply(M, x, xi, xnz, yscat, yi);

  //compress y
  spasm_compress_vector(yscat, yi, nz, y);

  //free workspace
  free(x);
  free(xi);
  free(yscat);

  return nz;
}
