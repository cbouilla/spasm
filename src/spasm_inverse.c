#include <assert.h>
#include "spasm.h"



/*
 * Given M, L, and k, compute y = x * M, 
 * where x is such as x * L = I[k].
 * the return value nz is the number of non-zero entries in y.
 */

int spasm_inverse_and_product(const spasm *L, const spasm *M, int k, spasm_GFp *y, int *yi, const int *pinv) {
  spasm * I;
  spasm_GFp *x;
  int i, p, top, nz, *xi, *xxi, Ln, Mn, prime, Mprime, xnz;

  /* check inputs */
  if(L == NULL || M == NULL) {
    return 0;
  }

  Ln = L->n;
  Mn = M->n;

  assert(Ln == Mn );

  prime = L->prime;
  Mprime = M->prime;
  assert(prime == Mprime);

  /* First solve x * L = ek */
  // get workspace.
  x = malloc(Mn * sizeof(spasm_GFp));
  xi = malloc(3*Mn * sizeof(int));
  spasm_vector_zero(xi, 3*Mn);
  spasm_vector_zero(x, Mn);

  I = spasm_identity(Mn, prime); // <--- Identity matrix

  top = spasm_sparse_backward_solve(L, I, k, xi, x, pinv,0);

  spasm_csr_free(I); // <--- free extra-workspace

  /* find x pattern, xxi. */
  xnz = Mn - top;
  xxi = malloc(xnz * sizeof(int));
  i = 0;
  for(p = top; p < Mn; p++) {
    xxi[i] = xi[p];
    i++; 
  }
  free(xi); // <--- free extra-workspace
 
  /* Calculate product */
  nz =  spasm_sparse_vector_matrix_prod(M, x, xxi, xnz, y, yi);

  /* free workspace */
  free(x);
  free(xxi);

  return nz;

}


/*
 * Given M, L, A and k, compute y = x * M, 
 * where x is such as x * L = A[k].
 * the return value nz is the number of non-zero entries in y.
 */
int spasm_solve_and_product(const spasm *L, const spasm *M, const spasm *A, int k, spasm_GFp *y, int *yi, const int *pinv) {
  spasm_GFp *x;
  int i, p, top, nz, *xi, *xxi, Ln, Mn, prime, xnz, Am;

  /* check inputs */
  if(L == NULL || M == NULL) {
    return 0;
  }

  Ln = L->n;
  Mn = M->n;
  Am = A->m;

  assert(Ln == Mn );
  assert(Am == Mn);

  prime = L->prime;
  assert(prime == M->prime);
  assert(prime == A->prime);

  /* First solve x * L = ek */
  // get workspace.
  x = malloc(Mn * sizeof(spasm_GFp));
  xi = malloc(3*Mn * sizeof(int));
  spasm_vector_zero(xi, 3*Mn);
  spasm_vector_zero(x, Mn);

  top = spasm_sparse_backward_solve(L, A, k, xi, x, pinv,0);

  /* find x pattern, xxi. */
  xnz = Mn - top;
  xxi = malloc(xnz * sizeof(int));
  i = 0;
  for(p = top; p < Mn; p++) {
    xxi[i] = xi[p];
    i++; 
  }
  free(xi); // <--- free extra-workspace
 
  /* Calculate product */
  nz =  spasm_sparse_vector_matrix_prod(M, x, xxi, xnz, y, yi);

  /* free workspace */
  free(x);
  free(xxi);

  return nz;

}
