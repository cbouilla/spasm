#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

/*
 * compute a (P)L(QU) decomposition.
 *
 * r = min(n, m) is an upper-bound on the rank of A
 *
 * L n * r
 * U is r * m
 *
 * pinv[j] = i if the pivot on column j is on row i. -1 if no pivot (yet) found
 * on column j.
 *
 */
spasm_lu *spasm_LU(const spasm * A) {
    spasm *L, *U;
    spasm_lu *N;
    spasm_GFp *Lx, *Ux, *x;
    int *Lp, *Lj, *Up, *Uj, *pinv, *xi;
    int n, m, r, ipiv, i, j, top, p, lnz, unz, prime, defficiency;

    /* check inputs */
    assert(A != NULL);

    n = A->n;
    m = A->m;
    r = spasm_min(n, m);
    prime = A->prime;
    defficiency = 0;

    /* educated guess of the size of L,U */
    lnz = 4 * spasm_nnz(A) + n;
    unz = 4 * spasm_nnz(A) + n;

    /* get GFp workspace */
    x = spasm_malloc(m * sizeof(spasm_GFp));

    /* get int workspace */
    xi = spasm_malloc(3 * m * sizeof(int));

    /* allocate result */
    N = spasm_malloc(sizeof(spasm_lu));

    /* allocate result L */
    N->L = L = spasm_csr_alloc(n, r, lnz, prime, true);

    /* allocate result U */
    N->U = U = spasm_csr_alloc(r, m, unz, prime, true);

    /* allocate result pinv */
    N->pinv = pinv = spasm_malloc(m * sizeof(int));

    Lp = L->p;
    Up = U->p;

    /* clear workspace */
    for (i = 0; i < m; i++) {
        x[i] = 0;
    }

    /* no rows pivotal yet */
    for (i = 0; i < m; i++) {
        pinv[i] = -1;
    }

    /* no rows of U yet */
    for (i = 0; i <= r; i++) {
        Up[i] = 0;
    }
    lnz = unz = 0;

    /* compute L[i] and U[i] */
    for (i = 0; i < n; i++) {
        /* --- Triangular solve: x * U = A[i] ---------------------------------------- */
        Lp[i] = lnz;            /* L[i] starts here */
	Up[i - defficiency] = unz;            /* U[i] starts here */

        /* not enough room in L/U ? realloc twice the size */
        if (lnz + m > L->nzmax) {
            spasm_csr_realloc(L, 2 * L->nzmax + m);
        }
        if (unz + m > U->nzmax) {
            spasm_csr_realloc(U, 2 * U->nzmax + m);
        }
        Lj = L->j;
        Lx = L->x;
        Uj = U->j;
        Ux = U->x;

        top = spasm_sparse_forwardsolve(U, A, i, xi, x, pinv);

        /* Find pivot and dispatch coeffs into L and U --------------- */
        ipiv = -1;
        /* index of best pivot so far.*/

	for (p = top; p < m; p++) {
            /* x[j] is (generically) nonzero */
            j = xi[p];

            /* if x[j] == 0 (numerical cancelation), we just ignore it */
            if (x[j] == 0) {
                continue;
            }

            if (pinv[j] < 0) {
                /* column j is not yet pivotal ? */

                /* have found the pivot on row i yet ? */
                if (ipiv == -1 || j < ipiv) {
                    ipiv = j;
                }
            } else {
                /* column j is pivotal */
                /* x[j] is the entry L[i, pinv[j] ] */
                Lj[lnz] = pinv[j];
                Lx[lnz] = x[j];
                lnz++;
            }
        }

        /* pivot found */
        if (ipiv != -1) {

	  /* L[i,i] <--- 1 */
	  Lj[lnz] = i - defficiency;
	  Lx[lnz] = 1;
	  lnz++;
	  pinv[ ipiv ] = i - defficiency;

	  /* pivot must be the first entry in U[i] */
	  Uj[unz] = ipiv;
	  Ux[unz] = x[ ipiv ];
	  unz++;

	  /* send remaining non-pivot coefficients into U */
	  for (p = top; p < m; p++) {
	    j = xi[p];

	    if (pinv[j] < 0) {
	      Uj[unz] = j;
	      Ux[unz] = x[j];
	      unz++;
	    }
	  }
	} else {
	  defficiency++;
	}
    }

    /* --- Finalize L and U ------------------------------------------------- */
    Lp[n] = lnz;
    Up[n - defficiency] = unz;
    spasm_csr_resize(U, n - defficiency, m);
    spasm_csr_resize(L, n, n - defficiency);

    /* remove extra space from L and U */
    spasm_csr_realloc(L, -1);
    spasm_csr_realloc(U, -1);
    free(x);
    free(xi);

    return N;
}

void spasm_free_LU(spasm_lu *X) {
  assert( X != NULL );
  spasm_csr_free(X->L);
  spasm_csr_free(X->U);
  free(X->pinv);
  free(X);
}
