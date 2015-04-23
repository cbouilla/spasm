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
 * Given M, L, and k, calculate y = x * M, 
 * where x is such as x * L = ek.
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
  if(Ln != Mn) {
    printf("Error : incompatible matrix size, %d, %d\n", Ln, Mn);
    return 0;
  }

  prime = L->prime;
  Mprime = M->prime;
  if(prime != Mprime) {
    printf("Error : incompatible matrix field, GF%d, GF%d\n", prime, Mprime);
    return 0;
  }

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


/* Given M, L, and two integers "from" and "to", calculate the product
 * Y = X * M, where row i of X, x_i, is such as x_i * L = e_(i+from).
 *
 * Y is given as a "spasm" matrix.
 */

void linvxm(const spasm *L, const spasm *M, int from, int to, spasm *Y, int *pinv) {
  int i, p, val, k, j, *yi, nz, *Yp, *Yj;
  spasm_GFp *y, m, *Yx;

  m = M->m;
  Yp = Y->p;
  Yj = Y->j;
  Yx = Y->x;

  // get workspace.
  y = malloc(m * sizeof(spasm_GFp));
  yi = malloc(m * sizeof(int));
  spasm_vector_zero(y, m);
  spasm_vector_zero(yi, m);

  // variables initialisation.
  Yp[0] = 0;
  p = 1;
  val = 0;
 
  for(i = from; i < to; i++) {

    // find the ith row of Y :
    nz = spasm_inverse_and_product(L, M, i, y, yi, pinv); // <--- nz : number of non zero-entries on row i.
    Yp[p] = Yp[p-1] + nz; // <--- update row pointers tab.
    p++;

    for(k = 0; k < nz; k++) {
      j = yi[k];     
      Yj[val] = j; // <--- add new column indice 
      Yx[val] = y[j]; // <--- add new numerical value
      val++;
    }

    // re-initialize y and yi.
    spasm_vector_zero(y, m); 
    spasm_vector_zero(yi, m);
  }

  // free workspace
  free(y);
  free(yi);
 
}
