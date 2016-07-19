/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_lu.c  */
#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

#ifdef SPASM_TIMING
uint64_t data_shuffling = 0;
#endif

#define DEBUG

/*
 * /!\ the ``row_permutation'' argument is NOT embedded in P. This means that
 * : a) L is ***really*** lower-(triangular/trapezoidal) 
 * b) PLUQ = row_permutation*A
 *
 * TODO : change this.
 */
spasm_lu *spasm_PLUQ(const spasm * A, const int *row_permutation, int keep_L) {
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
  for (j = 0; j < m; j++) {
    if (qinv[j] == -1) {
      qinv[j] = m - k;
      k++;
    }
  }

  /*
   * permute the columns of U in place. U becomes really upper-trapezoidal.
   */
  for (i = 0; i < r; i++) {
    for (px = Up[i]; px < Up[i + 1]; px++) {
      Uj[px] = qinv[Uj[px]];
    }
  }

  if (keep_L) {
    /*
     * permute the rows of L (not in place). L becomes really
     * lower-trapezoidal.
     */
    LL = spasm_permute(L, N->p, SPASM_IDENTITY_PERMUTATION, SPASM_WITH_NUMERICAL_VALUES);
    N->L = LL;
    spasm_csr_free(L);
  }
  return N;
}

/**
 *   Computes a random linear combination of A[k:].
 *   returns TRUE iff it belongs to the row-space of U.
 *   This means that with proba >= 1-1/p, all pivots have been found.
 */
int spasm_early_abort(const spasm * A, const int *row_permutation, int k, const spasm * U, int nu) {
  int *Aj, *Ap, *Uj, *Up;
  int i, j, inew, n, m, ok;
  spasm_GFp prime, *y, *Ax, *Ux;

  n = A->n;
  m = A->m;
  prime = A->prime;
  Aj = A->j;
  Ap = A->p;
  Ax = A->x;
  Uj = U->j;
  Up = U->p;
  Ux = U->x;

  y = spasm_malloc(m * sizeof(spasm_GFp));

  /* build a random (dense) linear combinations of the row of A, store in y */
  for (j = 0; j < m; j++) {
    y[j] = 0;
  }
  for (i = k; i < n; i++) {
    inew = (row_permutation != NULL) ? row_permutation[i] : i;
    spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, y, prime);
  }

  /* eliminate the coefficients of y below pivots of U */
  for (i = 0; i < nu; i++) {
    j = Uj[Up[i]];
    const spasm_GFp diagonal_entry = Ux[Up[i]];
    if (y[j] == 0) {
      continue;
    }
    const spasm_GFp d = (y[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
    spasm_scatter(Uj, Ux, Up[i], Up[i + 1], prime - d, y, prime);
  }

  /* if y != 0, then y does not belong to the row space of U */
  ok = 1;
  for (j = 0; j < m; j++) {
    if (y[j] != 0) {
      ok = 0;
      break;
    }
  }
  free(y);
  return ok;
}

/*
 * compute a (somewhat) LU decomposition using the GPLU algorithm.
 * 
 * r = min(n, m) is an upper-bound on the rank of A
 * 
 * L n * r U is r * m
 * 
L*U == row_permutation*A
 * 
 * qinv[j] = i if the pivot on column j is on row i. -1 if no pivot (yet) found
 * on column j.
 * 
 */
spasm_lu *spasm_LU(const spasm * A, const int *row_permutation, int keep_L) {
  spasm *L, *U;
  spasm_lu *N;
  spasm_GFp *Lx, *Ux, *x;
  int *Lp, *Lj, *Up, *Uj, *p, *qinv, *xj;
  int n, m, r, jpiv, i, inew, j, top, px, lnz, unz, old_unz, prime, defficiency,
      verbose_step;
  int rows_since_last_pivot, early_abort_done;

#ifdef SPASM_TIMING
  uint64_t start;
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
  xj = spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(xj, 3 * m);

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

  for (i = 0; i < n; i++) {
    /* no rows exchange yet */
    p[i] = i;
  }

  /* no rows of U yet */
  for (i = 0; i <= r; i++) {
    Up[i] = 0;
  }
  old_unz = lnz = unz = 0;

  /* initialize early abort */
  rows_since_last_pivot = 0;
  early_abort_done = 0;

  /* --- Main loop : compute L[i] and U[i] ------------------- */
  for (i = 0; i < n; i++) {
    if (!keep_L && i - defficiency == r) {
      fprintf(stderr, "\n[LU] full rank reached ; early abort\n");
      break;
    }
    if (!keep_L && !early_abort_done && rows_since_last_pivot > 10 && (rows_since_last_pivot > (n / 100))) {
      fprintf(stderr, "\n[LU] testing for early abort...");
      fflush(stderr);
      if (spasm_early_abort(A, row_permutation, i + 1, U, i - defficiency)) {
        fprintf(stderr, "SUCCESS\n");
        break;
      } else {
        fprintf(stderr, "FAILED\n");
        fflush(stderr);
      }
      early_abort_done = 1;
    }
    /*
     * --- Triangular solve: x * U = A[i]
     * ----------------------------------------
     */
    if (keep_L) {
      Lp[i] = lnz;              /* L[i] starts here */
    }
    Up[i - defficiency] = unz;  /* U[i] starts here */

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
    top = spasm_sparse_forward_solve(U, A, inew, xj, x, qinv);


    /*
     * --- Find pivot and dispatch coeffs into L and U
     * --------------------------
     */
#ifdef SPASM_TIMING
    start = spasm_ticks();
#endif
    jpiv = -1; /* column index of best pivot so far. */
    

    for (px = top; px < m; px++) {
      /* x[j] is (generically) nonzero */
      j = xj[px];

      /* if x[j] == 0 (numerical cancelation), we just ignore it */
      if (x[j] == 0) {
        continue;
      }
      if (qinv[j] < 0) {
        /* column j is not yet pivotal ? */

        /* have found the pivot on row i yet ? */
        if (jpiv == -1 || j < jpiv) {
          jpiv = j;
        }
      } else if (keep_L) {
        /* column j is pivotal */
        /* x[j] is the entry L[i, qinv[j] ] */
        Lj[lnz] = qinv[j];
        Lx[lnz] = x[j];
        lnz++;
      }
    }

    /* pivot found ? */
    if (jpiv != -1) {
      old_unz = unz;

      /* L[i,i] <--- x[jpiv]. Last entry of the row ! */
      if (keep_L) {
        Lj[lnz] = i - defficiency;
        Lx[lnz] = x[jpiv];
        lnz++;
      }
      qinv[jpiv] = i - defficiency;
      p[i - defficiency] = i;

      /* pivot must be the first entry in U[i] */
      Uj[unz] = jpiv;
      Ux[unz] = 1;
      unz++;

      /* send remaining non-pivot coefficients into U */
      spasm_GFp beta = spasm_GFp_inverse(x[jpiv], prime);
      for (px = top; px < m; px++) {
        j = xj[px];

        if (qinv[j] < 0) {
          Uj[unz] = j;
          Ux[unz] = (x[j] * beta) % prime;
          unz++;
        }
      }

      /* reset early abort */
      rows_since_last_pivot = 0;
      early_abort_done = 0;
    } else {
      defficiency++;
      p[n - defficiency] = i;
      rows_since_last_pivot++;
    }

#ifdef SPASM_TIMING
    data_shuffling += spasm_ticks() - start;
#endif

    if ((i % verbose_step) == 0) {
      fprintf(stderr, "\rLU : %d / %d [|L| = %d / |U| = %d] -- current density= (%.3f vs %.3f) --- rank >= %d", i, n, lnz, unz, 1.0 * (m - top) / (m), 1.0 * (unz - old_unz) / m, i - defficiency);
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
  free(xj);
  return N;
}


void spasm_free_LU(spasm_lu * X) {
  assert(X != NULL);
  spasm_csr_free(X->L);
  spasm_csr_free(X->U);
  free(X->qinv);
  free(X->p);
  free(X);
}

/* Computes the Schur complement, by eliminating the pivots located on
* rows p[0] ... p[n_pivots-1] of input matrix A. The pivots must be the
* entries on the lines. This returns a sparse representation of S.
*/
spasm *spasm_schur(const spasm * A, const int *p, const int n_pivots) {
  spasm *S;
  int *Sp, *Sj, Sn, Sm, m, n, snz, px, *xi, i, inew, top, j, *qinv, *q,
      verbose_step, *Ap, *Aj;
  spasm_GFp *Sx, *x;

  /* check inputs */
  assert(A != NULL);
  n = A->n;
  m = A->m;
  assert(n >= n_pivots);
  assert(m >= n_pivots);

  /* Get Workspace */
  Sn = n - n_pivots;
  Sm = m - n_pivots;
  snz = 4 * (Sn + Sm);          /* educated guess */
  S = spasm_csr_alloc(Sn, Sm, snz, A->prime, SPASM_WITH_NUMERICAL_VALUES);
  qinv = spasm_malloc(m * sizeof(int));
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(int));
  verbose_step = spasm_max(1, n / 1000);

  Ap = A->p;
  Aj = A->j;
  Sp = S->p;
  Sj = S->j;
  Sx = S->x;

  spasm_vector_zero(xi, 3 * m);

  /* build qinv from A [i.e. the "cheap pivots"] */
  for (j = 0; j < m; j++) {
    qinv[j] = -1;
  }
  for (i = 0; i < n_pivots; i++) {
    inew = p[i];
    j = Aj[Ap[inew]];
    qinv[j] = inew;             /* (inew, j) is a pivot */
  }

  /*
   * q sends the non-pivotal columns of A to the columns of S. It is not the
   * inverse of qinv...
   */
  q = spasm_malloc(m * sizeof(int));
  i = 0;
  for (j = 0; j < m; j++) {
    if (qinv[j] < 0) {
      q[j] = i;
      i++;
    } else {
      q[j] = -1;
    }
  }

  /* ---- compute Schur complement ----- */
  snz = 0;                      /* non-zero in S */
  Sn = 0;                       /* rows in S */
  fprintf(stderr, "Starting Schur complement computation...\n");
  for (i = n_pivots; i < n; i++) {
    Sp[Sn] = snz;               /* S[i] starts here */

    /* not enough room in S ? realloc twice the size */
    if (snz + Sm > S->nzmax) {
      spasm_csr_realloc(S, 2 * S->nzmax + Sm);
      Sj = S->j;
      Sx = S->x;
    }
    /* triangular solve */
    inew = p[i];
    top = spasm_sparse_forward_solve(A, A, inew, xi, x, qinv);

    /* dispatch x in S */
    for (px = top; px < m; px++) {
      j = xi[px];

      if (x[j] == 0) {          /* if x[j] == 0 (numerical cancelation), we
                                 * just ignore it */
        continue;
      }
      /* send non-pivot coefficients into S */
      if (q[j] >= 0) {
        Sj[snz] = q[j];
        Sx[snz] = x[j];
        snz++;
      }
    }
    Sn++;

    if ((i % verbose_step) == 0) {
      fprintf(stderr, "\rSchur : %d / %d [S=%d * %d, %d NNZ] -- current density= (%.3f)", i, n, Sn, Sm, snz, 1.0 * snz / (1.0 * Sm * Sn));
      fflush(stderr);
    }
  }
  /* finalize S */
  fprintf(stderr, "\n");
  Sp[S->n] = snz;
  spasm_csr_realloc(S, -1);

  /* free extra workspace */
  free(qinv);
  free(q);
  free(x);
  free(xi);

  return S;
}


/** Samples R rows at random in the schur complement of A w.r.t. the pivots in p[0:n_pivots],
* and return the number that are non-zero (these rows of A are linearly independent from the pivots).
*
* le seul truc utile, c'est la densité.
*/
int spasm_schur_probe(const spasm * A, const int *p, const int n_pivots, const int R, double *density) {
  int m, n, px, *xi, i, inew, top, j, *qinv, *Ap, *Aj, k, nnz, tmp;
  spasm_GFp *x;

  /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Aj = A->j;

  /* Get Workspace */
  qinv = spasm_malloc(m * sizeof(int));
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(xi, 3 * m);

  spasm_vector_set(qinv, 0, m, -1);
  for (i = 0; i < n_pivots; i++) {
    inew = p[i];
    j = Aj[Ap[inew]];
    qinv[j] = inew;             /* (inew, j) is a pivot */
  }

  k = 0;
  nnz = 0;
  for (i = 0; i < R; i++) {
    /* pick a random row in S, check if non-zero */
    inew = p[n_pivots + (rand() % (n - n_pivots))];
    top = spasm_sparse_forward_solve(A, A, inew, xi, x, qinv);
    tmp = 0;
    for (px = top; px < m; px++) {
      j = xi[px];
      if (qinv[j] < 0 && x[j] != 0) {  /* non-zero entry in a non-pivotal column ==> row of S in non-zero */
        nnz++;
        tmp = 1;
      }
    }
    k += tmp;
  }
    
  /* free extra workspace */
  free(qinv);
  free(x);
  free(xi);

  *density = ((double) nnz) / (m - n_pivots) / R;
  return (int) (((double) k) / ((double) R) * (n - n_pivots));
}


/* eliminate everything in the (dense) vector x using the pivots found in A */
void spasm_eliminate_sparse_pivots(const spasm * A, const int npiv, const int *p, spasm_GFp *x) {
  int i, inew, j, prime, *Aj, *Ap, k;
  spasm_GFp *Ax;

  Aj = A->j;
  Ap = A->p;
  Ax = A->x;
  prime = A->prime;

  for(i = 0; i < npiv; i++) {
    inew = p[i];
    j = Aj[Ap[inew]];
    
    if (x[j] == 0) {
      continue;
    }

    /* computing this inverse here is bad. The pivots in A should be unitary */
    // const spasm_GFp diagonal_entry = Ax[ Ap[inew] ];
    // assert (diagonal_entry == 1);

    
    spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], prime - x[j], x, prime);
    assert(x[j] == 0);

    for(k = 0; k < i; k++) {
      j = Aj[Ap[inew]];
      // assert(x[k] == 0); /* it stays eliminated */
    }
  }
}




int spasm_schur_rank(const spasm * A, const int *p, const int npiv) {
  int Sn, Sm, m, n, i, inew, j, k, r, px, prime, new, step, nbad, ngood;
  int *qinv, *q, *Ap, *Aj;
  spasm_GFp *Ax, *x, *y, *Uq;

  n = A->n;
  m = A->m;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  prime = A->prime;

  /* Get Workspace */
  Sn = n - npiv;
  Sm = m - npiv;
  q = spasm_malloc(Sm * sizeof(int));
  qinv = spasm_malloc(m * sizeof(int));
  x = spasm_malloc(m * sizeof(spasm_GFp));
  y = spasm_malloc(Sm * sizeof(spasm_GFp));
  
  /* build qinv from A and p [this sucks. qinv should be given] */
  spasm_vector_set(qinv, 0, m, -1);
  for (i = 0; i < npiv; i++) {
    inew = p[i];
    j = Aj[Ap[inew]];
    qinv[j] = inew;             /* (inew, j) is a pivot */
  }

  /* q sends columns of S to non-pivotal columns of A */
  k = 0;
  for (j = 0; j < m; j++) {
    if (qinv[j] < 0) {
      q[k] = j;
      k++;
    }
  }

  /* make pivotal rows of A unitary */
  for(i = 0; i < npiv; i++) {
    inew = p[i];
    const spasm_GFp diagonal_entry = Ax[ Ap[inew] ];
    const spasm_GFp alpha = spasm_GFp_inverse(diagonal_entry, prime);
    for (px = Ap[inew]; px < Ap[inew+1]; px++) {
      Ax[px] *= alpha;
      Ax[px] = Ax[px] % prime;
    }
  }

  spasm_dense_lu *U = spasm_dense_LU_alloc(Sm, A->prime);
  Uq = spasm_malloc(Sm * sizeof(int));
  for (j = 0; j < Sm; j++) {
    Uq[j] = j;
  }

  /* ---- compute Schur complement ----- */
  fprintf(stderr, "Starting Schur complement computation...\n");
  r = 0;
  step = 1;
  k = 1;
  nbad = 0;
  ngood = 0;
  int64 a=0, b=0, c=0, start;

  while (1) {
    /* compute a random linear combination of [step] non-pivotal rows. not DRY */
    
    /* l'essentiel du temps se passe avec step == 1 */

    start = spasm_ticks();
    spasm_vector_zero(x, m);
    if (step >= n) { /* on prend tout ! */
      for(i = npiv; i < n; i++) {
        inew = p[i];
        spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, x, prime);
      }
    } else { /* on prend un sous-ensemble aléatoire */
      for(i = 0; i < step; i++) {
        inew = p[npiv + (rand() % (n - npiv))];
        spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], 1 + (rand() % (prime-1)), x, prime);
      }
    }
    a += spasm_ticks() - start;
    start = spasm_ticks();

    spasm_eliminate_sparse_pivots(A, npiv, p, x);
    
    b += spasm_ticks() - start;
    start = spasm_ticks();

    /* the solution is scattered in x. Copy it to a (contiguous) y */
    for (j = 0; j < Sm; j++) {
      y[j] = x[ q[ Uq[j] ] ];
    }

    new = spasm_dense_LU_process(U, y, Uq);

    c += spasm_ticks() - start;

    if (new) {
      r++;
      ngood++;
      if (ngood == 16) {
        if (step > 1) {
          nbad = 0;
        }
        step = spasm_max(step / 2, 1);
        ngood = 0;
      }
      // fprintf(stderr, "+");
    } else {
      nbad++;
      if (nbad == 3) {
        if (step > n) {
          break; // c'est fini
        }
        step *= 2;
        nbad = 0;
        ngood = 0;
      }
      // fprintf(stderr, ".");
    }
    
    if ((k % 1) == 0) {
      fprintf(stderr, "\rSchur : %d [step = %d] -- current rank = %d / bad = %d, good = %d | %"PRId64" %"PRId64" %"PRId64"          ", k, step, r, nbad, ngood, a/k,b/k,c/k);
      fflush(stderr);
    }
    k++;
  }
  fprintf(stderr, "\n");
  return r;
}
