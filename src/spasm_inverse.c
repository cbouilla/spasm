#include <assert.h>
#include "spasm.h"


/*
 * Calculate y = x * M, where x and M are sparse.
 *
 * The result is scattered in y, its pattern is given by yi. 
 * The return value nz is the number of non-zero entries in y.
 */
int spasm_sparse_vector_matrix_prod(const spasm *M, const spasm_GFp *x, const int *xi, int xnz, spasm_GFp *y, int *yi) {
  int p, i, j, k, m, nz, Mnz=spasm_nnz(M), prime, *Mp, *Mj, *w;
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
  w = spasm_calloc(m, sizeof(int));

  /* primo, trouver support du résultat */
  nz = 0;
  for (k = 0; k < xnz; k++) {
    i = xi[k];

    for (p = Mp[i]; p < Mp[i+1]; p++) {
      j = Mj[p];
    
      if (w[j] == 0) {
	w[j] = 1;
	yi[nz] = j;
	nz++;
      }
    }
  }

  // scatter y
  for (k = 0; k < xnz; k++) {
    i = xi[k];

    spasm_scatter(Mj, Mx, Mp[i], Mp[i+1], x[i], y, prime);
  }

  // free workspace
  free(w);
  return nz;
}



/*
 * Given M, L, and k, compute y = x * M, 
 * where x is such as x * L = ek. (x is the k-th row of L^(-1))
 * the return value nz is the number of non-zero entries in y.
 */

int spasm_inverse_and_product(const spasm *L, const spasm *M, int k, spasm_GFp *y, int *yi, int *pinv) {
  spasm * I;
  spasm_GFp *x;
  int i, p, top, nz, *xi, *xxi, Ln, Mn, prime, Mprime, xnz;

  /* check inputs */
  if(L == NULL || M == NULL) {
    return 0;
  }

  Ln = L->n;
  Mn = M->n;
  assert(Ln = Mn );

  prime = L->prime;
  Mprime = M->prime;
  assert(prime = Mprime);


  /* First solve x * L = ek */
  // get workspace.
  x = malloc(Ln * sizeof(spasm_GFp));
  xi = malloc(3*Ln * sizeof(int));
  spasm_vector_zero(xi, 3*Ln);
  spasm_vector_zero(x, Ln);

  I = spasm_identity(Ln, prime); // <--- Identity matrix

  top = spasm_sparse_backward_solve(L, I, k, xi, x, pinv);

  free(I); // <--- free extra-workspace

  /* find x pattern, xxi. */
  xnz = Ln - top;
  xxi = malloc(xnz * sizeof(int));
  i = 0;
  for(p = top; p < Ln; p++) {
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


