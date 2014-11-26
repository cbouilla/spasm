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
void spasm_dense_backsolve(const spasm * L, spasm_GFp * x) {
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
 */
void spasm_dense_forwardsolve(const spasm * U, spasm_GFp * x) {
  int i, n, *Up, *Uj, prime;
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
      /* check diagonal entry */
      const spasm_GFp diagonal_entry = Ux[ Up[i] ];
      assert( diagonal_entry != 0 );

      // axpy - inplace
      x[i] = (x[i] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      spasm_scatter(Uj, Ux, Up[i] + 1, Up[i + 1], prime - x[i], x, prime);
    }
}



/*************** Triangular solving with sparse RHS
 *
 * solve G.x=B[:,k], where G is either upper (lo=0) or lower (lo=1) triangular.
 *
 * when this function returns, the solution scattered in x, and its pattern
 * is given in xi[top : n]. top is the return value.
 *
 */
int spasm_sparse_triangular_solve(spasm * G, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv, int lo) {
  int j, J, p, q, px, top, n, prime, *Gp, *Gj, *Bp, *Bj;
    spasm_GFp *Gx, *Bx;

    assert(G != NULL);
    assert(B != NULL);
    assert(xi != NULL);
    assert(x != NULL);

    n = G->n;
    Gp = G->p;
    Gj = G->j;
    Gx = G->x;
    prime = G->prime;

    Bp = B->p;
    Bj = B->j;
    Bx = B->x;

    /* xi[top : n] = Reach(B(:,k)) */
    top = spasm_reach(G, B, k, xi, pinv);

    /* clear x */
    for (p = top; p < n; p++) {
      x[xi[p]] = 0;
    }

    /* scatter B into x */
    for (p = Bp[k]; p < Bp[k + 1]; p++) {
        x[ Bj[p] ] = Bx[p];
    }

    /* iterate over the (precomputed) pattern of x (=the solution) */
    for (px = top; px < n; px++) {
      /* x(j) is nonzero */
      j = xi[px];

      /* j maps to col J of G */
      J = (pinv != NULL) ? (pinv[j]) : j;

      if (J < 0) {
	/* column J is empty */
            continue;
      }

      /* get G[j,j] */
      const spasm_GFp diagonal_entry = Gx[lo ? (Gp[J]) : (Gp[J + 1] - 1)];
      x[j] = (x[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      p = lo ? (Gp[J] + 1) : (Gp[J]);         /* lo: L(j,j) 1st entry */
      q = lo ? (Gp[J + 1]) : (Gp[J + 1] - 1); /* up: U(j,j) last entry */

      // this is a scatter
      spasm_scatter(Gj, Gx, p, q, prime - x[j], x, prime);
      //      for (; p < q; p++) {
      //	x[Gj[p]] -= Gx[p] * x[j];   /* x(i) -= G(i,j) * x(j) */
      //      }
    }
    return top;
}
