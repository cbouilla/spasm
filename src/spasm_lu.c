#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

#ifdef SPASM_TIMING
uint64_t data_shuffling = 0;
#endif

#define DEBUG

/*
 * /!\ the ``row_permutation'' argument is NOT embedded in P. This means that :
 *  a) L is ***really*** lower-(triangular/trapezoidal)
 *  b) PLUQ = row_permutation*A
 */
spasm_lu * spasm_PLUQ(const spasm *A, const int *row_permutation, int keep_L) {
  int m, i, j, r, px, k;
  int *Up, *Uj, *qinv;
  spasm *U, *L, *LL;
  spasm_lu *N;

  m = A->m;
  N = spasm_LU(A, row_permutation, keep_L);
  L = N->L;
  U = N->U;
  r = U->n;
  Up = U->p;
  Uj = U->j;
  qinv = N->qinv;

  k = 1;
  for(j = 0; j < m; j++) {
    if (qinv[j] == -1) {
      qinv[j] = m - k;
      k++;
    }
  }

  /* permute the columns of U in place.
     U becomes really upper-trapezoidal. */
  for(i = 0; i < r; i++) {
    for(px = Up[i]; px < Up[i + 1]; px++) {
      Uj[px] = qinv[ Uj[px] ];
    }
  }

  if (keep_L) {
    /* permute the rows of L (not in place).
       L becomes really lower-trapezoidal. */
    LL = spasm_permute(L, N->p, SPASM_IDENTITY_PERMUTATION, SPASM_WITH_NUMERICAL_VALUES);
    N->L = LL;
    spasm_csr_free(L);
  }
  
  return N;
}

/*
 * compute a (somewhat) LU decomposition.
 *
 * r = min(n, m) is an upper-bound on the rank of A
 *
 * L n * r
 * U is r * m
 *
 * qinv[j] = i if the pivot on column j is on row i. -1 if no pivot (yet) found
 * on column j.
 *
 */
spasm_lu *spasm_LU(const spasm * A, const int *row_permutation, int keep_L) {
    spasm *L, *U;
    spasm_lu *N;
    spasm_GFp *Lx, *Ux, *x;
    int *Lp, *Lj, *Up, *Uj, *p, *qinv, *xi;
    int n, m, r, ipiv, i, inew, j, top, px, lnz, unz, old_unz, prime, defficiency, verbose_step;

#ifdef SPASM_TIMING
    uint64_t start;
#endif

#ifdef SPASM_COL_WEIGHT_PIVOT_SELECTION
    int piv_weight;
    int *col_weights, *Ap, *Aj;
#endif

    /* check inputs */
    assert(A != NULL);

    n = A->n;
    m = A->m;
    r = spasm_min(n, m);
    prime = A->prime;
    defficiency = 0;
    verbose_step = spasm_max(1, n / 1000);

    /* educated guess of the size of L,U */
    lnz = 4 * spasm_nnz(A) + n;
    unz = 4 * spasm_nnz(A) + n;

    /* get GFp workspace */
    x = spasm_malloc(m * sizeof(spasm_GFp));

    /* get int workspace */
    xi = spasm_malloc(3 * m * sizeof(int));
    spasm_vector_zero(xi, 3*m);

    /* allocate result */
    N = spasm_malloc(sizeof(spasm_lu));
    N->L = L = (keep_L) ? spasm_csr_alloc(n, r, lnz, prime, true) : NULL;
    N->U = U = spasm_csr_alloc(r, m, unz, prime, true);
    N->qinv = qinv = spasm_malloc(m * sizeof(int));
    N->p = p = spasm_malloc(n * sizeof(int));

    Lp = (keep_L) ? L->p : NULL;
    Up = U->p;

    for (i = 0; i < m; i++) {
      /* clear workspace */
      x[i] = 0;
      /* no rows pivotal yet */
      qinv[i] = -1;
    }

#ifdef SPASM_COL_WEIGHT_PIVOT_SELECTION
    col_weights = spasm_malloc(m * sizeof(int));
    Ap = A->p;
    Aj = A->j;
    for (i = 0; i < m; i++) {
      col_weights[i] = 0;
    }
    for (i = 0; i < n; i++) {
      for(px = Ap[i]; px < Ap[i + 1]; px++) {
	col_weights[ Aj[px] ]++;
      }
    }
#endif

    for (i = 0; i < n; i++) {
      /* no rows exchange yet */
      p[i] = i;
    }

    /* no rows of U yet */
    for (i = 0; i <= r; i++) {
        Up[i] = 0;
    }
    old_unz = lnz = unz = 0;

    /* --- Main loop : compute L[i] and U[i] ------------------- */
    for (i = 0; i < n; i++) {
      if (!keep_L && i - defficiency == r) {
	fprintf(stderr, "\n[LU] full rank reached ; early abort\n");
	break;
      }
 
        /* --- Triangular solve: x * U = A[i] ---------------------------------------- */
      if (keep_L) {
	Lp[i] = lnz;                          /* L[i] starts here */
      }
      Up[i - defficiency] = unz;            /* U[i] starts here */

        /* not enough room in L/U ? realloc twice the size */
        if (keep_L && lnz + m > L->nzmax) {
            spasm_csr_realloc(L, 2 * L->nzmax + m);
        }
        if (unz + m > U->nzmax) {
            spasm_csr_realloc(U, 2 * U->nzmax + m);
        }
        Lj = (keep_L) ? L->j : NULL;
        Lx = (keep_L) ? L->x : NULL;
        Uj = U->j;
        Ux = U->x;

	inew = (row_permutation != NULL) ? row_permutation[i] : i;
        top = spasm_sparse_forward_solve(U, A, inew, xi, x, qinv);


        /* --- Find pivot and dispatch coeffs into L and U -------------------------- */
#ifdef SPASM_TIMING
      start = spasm_ticks();
#endif
        ipiv = -1;
        /* index of best pivot so far.*/

	for (px = top; px < m; px++) {
            /* x[j] is (generically) nonzero */
            j = xi[px];

            /* if x[j] == 0 (numerical cancelation), we just ignore it */
            if (x[j] == 0) {
                continue;
            }

            if (qinv[j] < 0) {
                /* column j is not yet pivotal ? */

                /* have found the pivot on row i yet ? */
#ifdef SPASM_COL_WEIGHT_PIVOT_SELECTION
                if (ipiv == -1 || col_weights[j] < piv_weight ) {
                    ipiv = j;
		    piv_weight = col_weights[j];
                }
#else
                if (ipiv == -1 || j < ipiv) {
                    ipiv = j;
                }
#endif
            } else if (keep_L) {
                /* column j is pivotal */
                /* x[j] is the entry L[i, qinv[j] ] */
                Lj[lnz] = qinv[j];
                Lx[lnz] = x[j];
                lnz++;
            }
        }

        /* pivot found */
        if (ipiv != -1) {
	  old_unz = unz;
	  //	  printf("\n pivot found on row %d of A at column %d\n", inew, ipiv);

	  /* L[i,i] <--- 1. Last entry of the row ! */
	  if (keep_L) {
	    Lj[lnz] = i - defficiency;
	    Lx[lnz] = 1;
	    lnz++;
	  }

	  qinv[ ipiv ] = i - defficiency;
	  p[i - defficiency] = i;

	  /* pivot must be the first entry in U[i] */
	  Uj[unz] = ipiv;
	  Ux[unz] = x[ ipiv ];
	  unz++;

	  /* send remaining non-pivot coefficients into U */
	  for (px = top; px < m; px++) {
	    j = xi[px];

	    if (qinv[j] < 0) {
	      Uj[unz] = j;
	      Ux[unz] = x[j];
	      unz++;
	    }
	  }

	} else {
	  defficiency++;
	  p[n - defficiency] = i;
	}


#ifdef SPASM_COL_WEIGHT_PIVOT_SELECTION
	/* update remaining col weights */
	for(px = Ap[inew]; px < Ap[inew + 1]; px++) {
	  col_weights[ Aj[px] ]--;
	}
#endif


#ifdef SPASM_TIMING
      data_shuffling += spasm_ticks() - start;
#endif

      if ((i % verbose_step) == 0) {
	fprintf(stderr, "\rLU : %d / %d [|L| = %d / |U| = %d] -- current density= (%.3f vs %.3f) --- rank >= %d", i, n, lnz, unz, 1.0 * (m-top) / (m), 1.0 * (unz-old_unz) / m, i - defficiency);
	fflush(stderr);
      }

    }
    
    /* --- Finalize L and U ------------------------------------------------- */
    fprintf(stderr, "\n");

    /* remove extra space from L and U */
    Up[i - defficiency] = unz;
    spasm_csr_resize(U, i - defficiency, m);
    spasm_csr_realloc(U, -1);

    if (keep_L) {
      Lp[n] = lnz;
      spasm_csr_resize(L, n, n - defficiency);
      spasm_csr_realloc(L, -1);
    }

    free(x);
    free(xi);

#ifdef SPASM_COL_WEIGHT_PIVOT_SELECTION
    free(col_weights);
#endif

    
    return N;
}


void spasm_free_LU(spasm_lu *X) {
  assert( X != NULL );
  spasm_csr_free(X->L);
  spasm_csr_free(X->U);
  free(X->qinv);
  free(X->p);
  free(X);
}


/*
 * Dispatch x entries in U (and L), find pivot if there is one.
 *
 * xi : x pattern.
 * top : first index in xi.
 * unz : number of entries in U. lnz : number of entries in L.
 * deff : rank defficiency. Current row in U is i - def.
 * i : current row in L (also current row in U if def = 0)
 * p : rows permutation. qinv : inverse of columns permutation.
 * n : size of p.
 *
 * return value is 1 if a pivot has been found, and 0 otherwise.
 */
int spasm_find_pivot(int *xi, spasm_GFp *x, int top, spasm *U, spasm *L, int *unz_ptr, int *lnz_ptr, int i, int *deff_ptr, int *qinv, int *p, int n) {
  int found, unz, lnz, m, ipiv, j, px, deff;
  int *Uj, *Lj;
  spasm_GFp *Ux, *Lx;

  assert(U != NULL);

  m = U->m;

  Uj = U->j;
  Ux = U->x;
  Lj = (L != NULL) ? L->j : NULL;
  Lx = (L != NULL) ? L->x : NULL;

  unz = *unz_ptr;
  lnz = *lnz_ptr;
  deff = *deff_ptr;
 
  /* --- Find pivot and dispatch coeffs ----------------------------------- */
  ipiv = -1; //<--- no pivot found so far.

  //search pivot in genericaly nonzero entries of x.
  for(px = top; px < m; px++) {
    j = xi[px]; 
    //if x[j] = 0 (numerical cancelation) we ignore it
    if(x[j] == 0) continue;

    if(qinv[j] < 0) {
      // no pivot on column j yet.
      
      if(ipiv == -1 || j < ipiv) {
	ipiv = j; // <--- best pivot so far is on column j.
      }

    }
    else if(L != NULL) {
      // column j already pivotal.
      // dispatch x[j] in L[i_l, qinv[j]]
      Lj[lnz] = qinv[j];
      Lx[lnz] = x[j];
      lnz ++;
    }

  }

  // pivot found : 
  if(ipiv != -1) {
  
    found = 1;
    // Last entry on row i of L is 1.
    if(L != NULL) {
      Lj[lnz] = i - deff; 
      Lx[lnz] = 1;
      lnz++;
    } 

    // update permutation
    qinv[ipiv] = i - deff;
    p[i - deff] = i; 

    // add entries in U
    Uj[unz] = ipiv;
    Ux[unz] = x[ipiv]; // <--- add pivot first.
    unz++;

    for(px = top; px < m; px++) {
      // dispatch other entries in U.
      j = xi[px]; // <--- non zero entries

      //if(x[j] == 0) continue; //<-- if numerical cancelation, we ignore it.

      if(qinv[j] < 0) {
	// no pivot in column j yet
	Uj[unz] = j;
	Ux[unz] = x[j];
	unz++;
      }
    }

  }

  //no pivot found :
  else { 
    found = 0;
    deff++;
    p[n - deff] = i;
  }

  /* --- Clear workspace -------------------------------------------------- */
  spasm_vector_zero(xi, m);
  spasm_vector_zero(x, m);

  /* --- Return ----------------------------------------------------------- */
  *unz_ptr = unz; 
  *lnz_ptr = lnz;
  *deff_ptr = deff;

  return found;
}
