#include <assert.h>
#include "spasm.h"


/*
 * Calculate y = x * M, where x and M are sparse.
 *
 * The result is scattered in y, its pattern is given by yi. 
 * The return value ytop is the number of non-zero entries in y.
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
 * Given a matrix M and a lower_triangular matrix L, find y such as
 * y = x * M where x is solution of x * L = ek.
 *
 * the return value nz is the number of non-zero entries in y.
 */

int spasm_inverse_and_sparse_vect_prod(const spasm *L, const spasm *M, int k, spasm_GFp *y, int *yi, const int *pinv) {
  spasm *I;
  int Ln, prime, top, nz, *xi;
  spasm_GFp *x;

  Ln = L->n;
  prime = L->prime;

  /* First : find x such as x * L = ek */
  I = spasm_identity(n, prime); // <---- identity matrix
  // get workspace
  x = malloc(Ln * sizeof(spasm_GFp));
  xi = malloc(3*Ln * sizeof(int));
  spasm_vector_zero(xi, 3*n);
  spasm_vector_zero(x, n);

  // solve system, get top value.
  top = spasm_sparse_backward_solve(L, I, k, xi, x, pinv);

  spasm_csr_free(I); // free identity matrix




}
