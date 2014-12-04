#include <assert.h>
#include "spasm.h"
/* Solves x.A = b where A is not necessarily square;
 *
 * b has size m, solution has size n.
 *
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_PLUQ_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x) {
  spasm_GFp *y, *z;
    spasm_lu *PLUQ;
    spasm *L, *U;
    int n, m, ok ;

    /* check inputs */
    assert(A != NULL);
    assert(b != NULL);

    PLUQ = spasm_PLUQ(A);
    L = PLUQ->L;
    U = PLUQ->U;

    n = A->n;
    m = A->m;
    //    r = spasm_max(n, m);

    /* get workspace */
    y = spasm_malloc (m * sizeof(spasm_GFp));
    z = spasm_malloc (n * sizeof(spasm_GFp));

    /* y*Q = b */
    spasm_ipvec(PLUQ->qinv, b, y, m);

    /* z.U*Q = b  (if possible) */
    ok = spasm_dense_forward_solve(U, y, z, SPASM_IDENTITY_PERMUTATION);

    if (ok == SPASM_SUCCESS) {
      /* y.LUQ = b */
      spasm_dense_back_solve(L, z, y, SPASM_IDENTITY_PERMUTATION);

      /* x.PLUQ = b */
      spasm_ipvec(PLUQ->p, y, x, n);
    }

    free(y) ;
    free(z) ;
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
  spasm_GFp *y, *z, *foo;
    spasm_lu *LU;
    spasm *L, *U;
    int n, m, i, ok ;
    int *q;

    /* check inputs */
    assert(A != NULL);
    assert(b != NULL);

    LU = spasm_LU(A);
    L = LU->L;
    U = LU->U;

    n = A->n;
    m = A->m;
    //    r = spasm_max(n, m);

    /* get workspace */
    y = spasm_malloc(m * sizeof(spasm_GFp));
    z = spasm_malloc(U->n * sizeof(spasm_GFp));
    foo = spasm_malloc(m * sizeof(spasm_GFp));

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
      // defense
      for(i = 0; i < m; i++) {
	foo[i] = 0;
      }
      spasm_gaxpy(U, z, foo);
      for(i = 0; i < m; i++) {
	assert(foo[i] == b[i]);
      }
      // fin

      /* y.LU = b */
      spasm_dense_back_solve(L, z, x, LU->p);
      for(i = 0; i < m; i++) {
	foo[i] = 0;
      }
      spasm_gaxpy(A, x, foo);
      for(i = 0; i < m; i++) {
	assert(foo[i] == b[i]);
      }
    }

    free(y);
    free(z);
    free(q);
    spasm_free_LU(LU);
    return ok;
}
