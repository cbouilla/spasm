#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

/* compute a (P)L(QU) decomposition.
 *
 * L is always square of size n * n
 * U has the same size as A
 *
 * pinv[j] = i if the pivot on column j is on row i. -1 if no pivot (yet) found on column j.
 *
 */
spasm_lu *spasm_LU(const spasm * A) {
    spasm *L, *U;
    spasm_lu *N;
    spasm_GFp *Lx, *Ux, *x;
    int *Lp, *Lj, *Up, *Uj, *pinv, *xi;
    int n, m, ipiv, i, j, top, p, lnz, unz, prime;

    /* check inputs */
    assert(A != NULL);

    n = A->n;
    m = A->m;
    prime = A->prime;

    /* educated guess of the size of L,U */
    lnz = 4 * spasm_nnz(A) + n;
    unz = 4 * spasm_nnz(A) + n;

    /* get GFp workspace */
    x = spasm_malloc(n * sizeof(spasm_GFp));

    /* get int workspace */
    xi = spasm_malloc(2 * n * sizeof(int));

    /* allocate result */
    N = spasm_malloc( sizeof(spasm_lu) );

    /* allocate result L */
    N->L = L = spasm_csr_alloc(n, n, lnz, prime, true);

    /* allocate result U */
    N->U = U = spasm_csr_alloc(n, m, unz, prime, true);

    /* allocate result pinv */
    N->pinv = pinv = spasm_malloc(m * sizeof(int));

    Lp = L->p;
    Up = U->p;

    /* clear workspace */
    for (i = 0; i < n; i++) {
	x[i] = 0;
    }

    /* no rows pivotal yet */
    for (i = 0; i < m; i++) {
	pinv[i] = -1;
    }

    /* no rows of L yet */
    for (i = 0; i <= n; i++) {
	Lp[i] = 0;
    }
    lnz = unz = 0;

    /* compute L[i] and U[i] */
    for (i = 0; i < n; i++) {

      /* --- Triangular solve --------------------------------------------- */
	Lp[i] = lnz;            /* L[i] starts here */
	Up[i] = unz;            /* U[i] starts here */

	/* not enough room in L/U ? realloc twice the size */
	if (lnz + n > L->nzmax) {
	  spasm_csr_realloc(L, 2 * L->nzmax + n);
	}
	if (unz + n > U->nzmax) {
	  spasm_csr_realloc(U, 2 * U->nzmax + n);
	}
	Lj = L->j;
	Lx = L->x;
	Uj = U->j;
	Ux = U->x;
	// permutation des lignes ?

	/* Solve x * U = A[i] */
	top = spasm_sparse_forwardsolve(U, A, i, xi, x, pinv);

	/* --- Find pivot and dispatch coeffs into L and U ------------------------ */
	ipiv = -1; // index of best pivot so far.

	for (p = top; p < n; p++) {
	  /* x[j] is (generically) nonzero */
	  j = xi[p];

	  /* if x[j] == 0 (numerical cancelation), we just ignore it */
	  if (x[j] == 0) {
	    continue;
	  }

	    if (pinv[j] < 0) {
	      /* column j is not yet pivotal ? */
	      /* send coefficient into U[i,j] */
	      Uj[unz] = j;
	      Ux[unz] = x[j];
	      unz++;

	      /* have found the pivot on row i yet ? */
	      if (ipiv == -1) {
		ipiv = j;
		pinv[j] = i;
	      }
	    } else {
	      /* column j is pivotal */
	      /* x[j] is the entry L[i, pinv[j] ] */
		Lj[lnz] = pinv[j];
		Lx[lnz] = x[j];
		lnz++;
	    }

	    /* clean x for the next iteration */
	    x[j] = 0;
	}

	/* set L[i, i] = 1 if a pivot has been found */
	if (ipiv != -1 ) {
	  Lj[lnz] = i;
	  Lx[lnz] = 1;
	  lnz++;
	}
    }

    /* --- Finalize L and U ------------------------------------------------- */
    Lp[n] = lnz;
    Up[n] = unz;

    /* remove extra space from L and U */
    spasm_csr_realloc(L, 0);
    spasm_csr_realloc(U, 0);
    free(x);
    free(xi);

    return N;
}
