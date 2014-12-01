#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

void print_vec_(int *Aj, spasm_GFp * Ax, int p, int q, int n) {
    spasm_GFp x[n];
    int i, k;

    for (i = 0; i < n; i++) {
        x[i] = 0;
    }

    for (k = p; k < q; k++) {
        x[Aj[k]] = Ax[k];
    }

    printf("[ ");
    for (i = 0; i < n; i++) {
        printf("%d ", x[i]);
    }
    printf("]");
}

void print_vec(const spasm * A, int i) {
    print_vec_(A->j, A->x, A->p[i], A->p[i + 1], A->m);
}

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

    /*    printf("[DEBUG LU] input matrix is %d x %d modulo %d\n", n, m, prime);
    printf("[DEBUG LU] L will be %d x %d and U will be %d x %d\n", n, r, r, m);
    */
    
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
      /*      printf("[DEBUG LU] i = %d, rank defficiency = %d\n", i, defficiency);
      printf("A[i] = ");
        print_vec(A, i);
        printf("\n");
      */

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

	/*	for (p = 0; p < i - defficiency; p++) {
            printf("U[%d] = ", p);
            print_vec(U, p);
            printf("\n");
        }
	printf("-----------------\n");
	for (p = 0; p < i; p++) {
            printf("L[%d] = ", p);
            print_vec(L, p);
            printf("\n");
	    }*/

        /*  Solve  ------------- */

       	/*printf("[DEBUG LU] pinv = [");
        for (int r = 0; r < m; r++) {
	  printf("%d ", pinv[r]);
	}
	printf("]\n");
	*/

	/* bien nécessaire ? --> pour l'affichage seulement
	for(p = 0; p < m; p++) {
	  x[p] = 0;
	  }*/
        top = spasm_sparse_forwardsolve(U, A, i, xi, x, pinv);

	/* check resolution
	spasm_GFp y[m];
	for(p = 0; p < m; p++) {
	  y[p] = 0;
	}
	for(p = top; p < m; p++) {
	  j = xi[p];
	  if (pinv[j] == -1) {
	    y[j] = (y[j] + x[j]) % A->prime;
	  } else {
	    spasm_scatter(Uj, Ux, Up[ pinv[j] ], Up[pinv[j] + 1], x[j], y, prime);
	  }
	}
	spasm_scatter(A->j, A->x, A->p[i], A->p[i + 1], prime - 1, y, prime);

        printf("[DEBUG LU] x = [");
        for (int r = 0; r < m; r++) {
            printf("%d ", x[r]);
        }
        printf("]\n");
	//	printf("[DEBUG LU] y = [");
        for (int r = 0; r < m; r++) {
          //  printf("%d ", y[r]);
	  assert( y[r] == 0 );
        }
	//        printf("]\n");
	*/


        /*
         * --- Find pivot and dispatch coeffs into L and U
         * ------------------------
         */
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
	    /* clean x for the next iteration -- inutile fait par sparse_solve */
	    //	    x[j] = 0;
	  }
	} else {
	  defficiency++;
	}
	/*	printf("[DEBUG LU] pivot choisi colonne %d\n", ipiv);
        printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");
	*/
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
