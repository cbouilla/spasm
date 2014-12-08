#include <assert.h>
#include "spasm.h"
/* Solves x.A = b where A is not necessarily square;
 *
 * b has size m, solution has size n.
 *
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_PLUQ_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x) {
  spasm_GFp *u, *v, *w, *s;
    spasm_lu *PLUQ;
    spasm *L, *U;
    int n, m, r, ok;
    int *row_permutation;

    /* check inputs */
    assert(A != NULL);
    assert(b != NULL);

    row_permutation = spasm_row_sort(A);
    PLUQ = spasm_PLUQ(A, row_permutation);
    L = PLUQ->L;
    U = PLUQ->U;

    n = A->n;
    m = A->m;
    r = U->n;

    /* get workspace */
    u = spasm_malloc (m * sizeof(spasm_GFp));
    v = spasm_malloc (r * sizeof(spasm_GFp));
    w = spasm_malloc (n * sizeof(spasm_GFp));
    s = spasm_malloc (n * sizeof(spasm_GFp));

    /* u*Q = b */
    spasm_ipvec(PLUQ->qinv, b, u, m);

    /* v.U*Q = b  (if possible) */
    ok = spasm_dense_forward_solve(U, u, v, SPASM_IDENTITY_PERMUTATION);

    if (ok == SPASM_SUCCESS) {
      /* w.LUQ = b */
      spasm_dense_back_solve(L, v, w, SPASM_IDENTITY_PERMUTATION);

      /* x.PLUQ = b */
      spasm_ipvec(PLUQ->p, w, s, n);
      spasm_ipvec(row_permutation, s, x, n);
    }

    free(u);
    free(v);
    free(w);
    free(s);
    spasm_free_LU(PLUQ);
    return ok;
}

/* Solves x.A = b where A is not necessarily square;
 *
 * b has size m, solution has size n.
 *
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_LU_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x) {
  spasm_GFp *y, *z, *w;
  spasm_lu *LU;
  spasm *L, *U;
  int n, m, r, i, ok ;
  int *q, *row_permutation;

    /* check inputs */
    assert(A != NULL);
    assert(b != NULL);

    row_permutation = spasm_row_sort(A);
    LU = spasm_LU(A, row_permutation, SPASM_KEEP_L);
    L = LU->L;
    U = LU->U;

    n = A->n;
    m = A->m;
    r = U->n;

    /* get workspace */
    y = spasm_malloc(m * sizeof(spasm_GFp));
    z = spasm_malloc(r * sizeof(spasm_GFp));
    w = spasm_malloc(n * sizeof(spasm_GFp));
    q = spasm_malloc(m * sizeof(int));

    for(i = 0; i < m; i++) {
      if (LU->qinv[i] != -1) {
	q[ LU->qinv[i] ] = i;
      }
    }

    /* z.U = b  (if possible) */
    for(i = 0; i < m; i++) {
      y[i] = b[i];
    }
    ok = spasm_dense_forward_solve(U, y, z, q);

    if (ok == SPASM_SUCCESS) {

      /* y.LU = b */
      spasm_dense_back_solve(L, z, w, LU->p);
      spasm_ipvec(row_permutation, w, x, n);
    }

    free(y);
    free(z);
    free(w);
    free(q);
    free(row_permutation);
    spasm_free_LU(LU);
    return ok;
}
