#include <assert.h>
#include "spasm.h"

/*
 * Solving triangular systems, dense RHS
 */

// TODO : merge usolve and lsolve

/*
 * solve G.x=b where x and b are dense.  x=b on input, solution on output.
 *
 * G is lower-triangular (lo == 1) or upper-triangular (lo == 0)
 *
 * Assumption : the diagonal entry is always present, is always != 0.
 *
 * If G is lower-triangular, then it is the first one of each column.
 * More precisely, L[j,j] is Lx[ Lp[j] ]
 *
 * If G is upper-triangular, then it is the last one of each column.
 * More precisely, U[j,j] is Ux[ Up[j + 1] - 1 ]
 */
void spasm_triangular_solve(const spasm * G, spasm_GFp * x, int lo) {
  int p, q, j, n, from, to, step, *Gp, *Gi, prime;
    spasm_GFp *Gx;

    /* check inputs */
    assert(x != NULL);
    assert(spasm_is_csc(G));

    n = G->n;
    Gp = G->p;
    Gi = G->i;
    Gx = G->x;
    prime = G->prime;

    from = lo ? 0 : n-1;
    to   = lo ? n : -1;
    step = lo ? 1 : -1;

    for (j = from; j != to; j += step) {
      /* check diagonal entry */
      if (lo) {
	assert( Gi[ Gp[j] ] == j);
      } else {
	assert( Gi[ Gp[j + 1] - 1] == j );
      }
      const spasm_GFp diagonal_entry = lo ? Gx[ Gp[j] ] :  Gx[ Gp[j + 1] - 1 ];

      // axpy - inplace
      x[j] = (x[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;

      p = lo ? (Gp[j] + 1) : Gp[j];         /* lo: L(j,j) 1st entry */
      q = lo ? (Gp[j + 1]) : Gp[j + 1] - 1; /* up: U(j,j) last entry */

      spasm_scatter(Gi, Gx, p, q, prime - x[j], x, prime);
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
  int j, J, p, q, px, top, n, prime, *Gp, *Gi, *Bp, *Bi;
    spasm_GFp *Gx, *Bx;

    assert(spasm_is_csc(G));
    assert(spasm_is_csc(B));
    assert(xi != NULL);
    assert(x != NULL);

    n = G->n;
    Gp = G->p;
    Gi = G->i;
    Gx = G->x;
    prime = G->prime;

    Bp = B->p;
    Bi = B->i;
    Bx = B->x;

    /* xi[top : n] = Reach(B(:,k)) */
    top = spasm_reach(G, B, k, xi, pinv);

    /* clear x */
    for (p = top; p < n; p++) {
      x[xi[p]] = 0;
    }

    /* scatter B into x */
    for (p = Bp[k]; p < Bp[k + 1]; p++) {
        x[ Bi[p] ] = Bx[p];
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
      spasm_scatter(Gi, Gx, p, q, prime - x[j], x, prime);
      //      for (; p < q; p++) {
      //	x[Gi[p]] -= Gx[p] * x[j];   /* x(i) -= G(i,j) * x(j) */
      //      }
    }
    return top;
}
