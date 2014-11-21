#include <assert.h>
#include "spasm.h"

/*
 * Solving triangular systems, dense RHS
 */

// TODO : merge usolve and lsolve

/*
 * solve L.x=b where x and b are dense.  x=b on input, solution on output.
 *
 * Assumption : the diagonal entry is always present, is always != 0, and is the
 * first one of each column. More precisely, L[j,j] is Lx[ Lp[j] ]
 */
void spasm_lsolve(const spasm * L, spasm_GFp * x) {
    int p, j, n, *Lp, *Li, prime;
    spasm_GFp *Lx;

    /* check inputs */
    assert(x != NULL);
    assert(spasm_is_csc(L));

    n = L->n;
    Lp = L->p;
    Li = L->i;
    Lx = L->x;
    prime = L->prime;

    for (j = 0; j < n; j++) {
        /* check diagonal entry */
        assert(Li[Lp[j]] == j);
        const spasm_GFp diagonal_entry = Lx[ Lp[j] ];

        // axpy - inplace
	x[j] = (x[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
        spasm_scatter(Li, Lx, Lp[j] + 1, Lp[j + 1], prime - x[j], x, prime);
    }
}

/*
 * solve U.x=b where x and b are dense.  x=b on input, solution on output.
 *
 * Assumption : the diagonal entry is always present, is always != 0, and is the
 * last one of each column. More precisely, U[j,j] is Ux[ Up[j+1] - 1 ]
 */
void spasm_usolve(const spasm * U, spasm_GFp * x) {
    int p, j, n, *Up, *Ui, prime;
    spasm_GFp *Ux;

    /* check inputs */
    assert(x != NULL);
    assert(spasm_is_csc(U));

    n = U->n;
    Up = U->p;
    Ui = U->i;
    Ux = U->x;
    prime = U->prime;

    for (j = n - 1; j >= 0; j--) {
        /* check diagonal entry */
        assert(Ui[Up[j + 1] - 1] == j);
        const spasm_GFp diagonal_entry = Ux[Up[j + 1] - 1];

        //axpy - inplace
	x[j] = (x[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
        spasm_scatter(Ui, Ux, Up[j], Up[j + 1] - 1, prime - x[j], x, prime);
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
int spasm_spsolve(spasm * G, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv, int lo) {
    int j, J, p, q, px, top, n, *Gp, *Gi, *Bp, *Bi;
    spasm_GFp *Gx, *Bx;

    assert(spasm_is_csc(G));
    assert(spasm_is_csc(B));
    assert(xi != NULL);
    assert(x != NULL);

    n = G->n;
    Gp = G->p;
    Gi = G->i;
    Gx = G->x;

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
      spasm_scatter(Gi, Gx, p, q, x[j], x, prime);
      //      for (; p < q; p++) {
      //	x[Gi[p]] -= Gx[p] * x[j];   /* x(i) -= G(i,j) * x(j) */
      //      }
    }
    return top;
}
