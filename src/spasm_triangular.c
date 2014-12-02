#include <assert.h>
#include "spasm.h"

/*
 * Solving triangular systems, dense RHS
 */


/* dense backwards substitution solver. Solve x . L = b where x and b are dense.
 *
 * x=b on input, solution on output.
 *
 * L is assumed to be lower-triangular, with non-zero diagonal.
 *
 * The diagonal entry is the **last** of each column.
 * More precisely, L[j,j] is Lx[ Lp[j+1] - 1 ]
 */
void spasm_dense_back_solve(const spasm * L, spasm_GFp * x) {
  int i, n, *Lp, *Lj, prime;
    spasm_GFp *Lx;

    /* check inputs */
    assert(x != NULL);
    assert(L != NULL);

    n = L->n;
    Lp = L->p;
    Lj = L->j;
    Lx = L->x;
    prime = L->prime;

    for (i = n - 1; i >= 0; i--) {

      /* check diagonal entry */
      const spasm_GFp diagonal_entry = Lx[ Lp[i + 1] - 1 ];
      assert( diagonal_entry != 0 );

      // axpy - inplace
      x[i] = (x[i] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      spasm_scatter(Lj, Lx, Lp[i], Lp[i + 1] - 1, prime - x[i], x, prime);
    }
}

/* dense forwards substitution solver. Solve x . U = b where x and b are dense.
 *
 * x=b on input, solution on output.
 *
 * U is upper-triangular
 *
 * Assumption : the diagonal entry is always present, is always != 0.
 *
 * The diagonal entry is the first one of each row.
 * More precisely, U[i,i] is Ux[ Up[i] ]
 *
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_dense_forward_solve(const spasm * U, spasm_GFp * x, int* qinv) {
  int i, j, n, *Up, *Uj, prime;
    spasm_GFp *Ux;

    /* check inputs */
    assert(x != NULL);
    assert(U != NULL);

    n = U->n;
    Up = U->p;
    Uj = U->j;
    Ux = U->x;
    prime = U->prime;

    for (i = 0; i < n; i++) {
      if (x[i] != 0) {
	/* get pivot column */
	j = (qinv != NULL) ? qinv[i] : i;
	if (j < 0) {
	  return SPASM_NO_SOLUTION;
	}

	/* check diagonal entry */
	const spasm_GFp diagonal_entry = Ux[ Up[i] ];
	assert( diagonal_entry != 0 );

	// axpy - inplace
	x[i] = (x[i] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
	spasm_scatter(Uj, Ux, Up[i] + 1, Up[i + 1], prime - x[i], x, prime);
      }
    }
    return SPASM_SUCCESS;
}



/*************** Triangular solving with sparse RHS
 *
 * solve x * U = B[k], where U is (permuted) upper triangular.
 *
 * x has size m (number of columns of U, paradoxically).
 *
 * when this function returns, the solution scattered in x, and its pattern
 * is given in xi[top : m].
 *
 * top is the return value.
 *
 */
int spasm_sparse_forward_solve(spasm * U, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv) {
  int i, I, p, px, top, n, m, prime, *Up, *Uj, *Bp, *Bj;
    spasm_GFp *Ux, *Bx;

    assert(U != NULL);
    assert(B != NULL);
    assert(xi != NULL);
    assert(x != NULL);

    n = U->n;
    m = U->m;
    Up = U->p;
    Uj = U->j;
    Ux = U->x;
    prime = U->prime;

    Bp = B->p;
    Bj = B->j;
    Bx = B->x;

    /* xi[top : n] = Reach( U, B[k] ) */
    top = spasm_reach(U, B, k, xi, pinv);

    /* clear x */
    for (p = top; p < m; p++) {
      x[ xi[p] ] = 0;
    }

    /* scatter B[k] into x */
    for (p = Bp[k]; p < Bp[k + 1]; p++) {
        x[ Bj[p] ] = Bx[p];
    }

    /* iterate over the (precomputed) pattern of x (= the solution) */
    for (px = top; px < m; px++) {
      /* x[i] is nonzero */
      i = xi[px];

      /* i maps to row I of U */
      I = (pinv != NULL) ? (pinv[i]) : i;

      if (I < 0) {
	/* row I is empty */
            continue;
      }

      /* get U[i,i] */
      const spasm_GFp diagonal_entry = Ux[ Up[I] ];
      assert( diagonal_entry != 0 );
      // axpy-in-place
      x[i] = (x[i] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      spasm_scatter(Uj, Ux, Up[I] + 1, Up[I + 1], prime - x[i], x, prime);
    }
    return top;
}
