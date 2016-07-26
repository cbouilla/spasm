/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_lu.c  */
#include <assert.h>
#include "spasm.h"

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
*/
double spasm_schur_probe_density(const spasm * A, const int *p, const int *qinv, const int npiv, const int R) {
  int m, n, px, *xj, i, inew, top, j, *Ap, *Aj, nnz;
  spasm_GFp *x;

  /* check inputs */
  m = A->m;
  n = A->n;
  Ap = A->p;
  Aj = A->j;

  /* Get Workspace */
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xj = spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(xj, 3 * m);

  nnz = 0;
  for (i = 0; i < R; i++) {
    /* pick a random row in S, check if non-zero */
    inew = p[npiv + (rand() % (n - npiv))];
    top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);
    for (px = top; px < m; px++) {
      j = xj[px];
      if (qinv[j] < 0 && x[j] != 0) {  /* non-zero entry in a non-pivotal column ==> row of S in non-zero */
        nnz++;
      }
    }
  }
    
  /* free extra workspace */
  free(x);
  free(xj);

  return ((double) nnz) / (m - npiv) / R;
}


int spasm_schur_rank(const spasm * A, const int *p, const int npiv) {
  int Sn, Sm, m, n, i, inew, j, k, r, px, prime, new, step, nbad, ngood;
  int *qinv, *q, *Ap, *Aj;
  double wtime_start;
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
  fprintf(stderr, "rank of dense schur complement...\n");
  wtime_start = spasm_wtime();
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
  fprintf(stderr, "\n[schur/rank] Time: %.1fs\n", spasm_wtime() - wtime_start);
  return r;
}
