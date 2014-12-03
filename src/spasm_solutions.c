#include <assert.h>
#include "spasm.h"
/* Solves x.A = b where A is not necessarily square;

 *
 * b has size m, solution has size n.
 *
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_PLUQ_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x) {
  spasm_GFp *y;
    spasm_lu *PLUQ;
    spasm *L, *U;
    int i, n, m, r, ok ;

    /* check inputs */
    assert(A != NULL);
    assert(b != NULL);

    PLUQ = spasm_PLUQ(A);
    L = PLUQ->L;
    U = PLUQ->U;

    n = A->n;
    m = A->m;
    r = spasm_max(n, m);

    /* get workspace */
    y = spasm_malloc (r * sizeof(spasm_GFp));

    /* y*Q = b */
    spasm_ipvec(PLUQ->qinv, b, y, m);

    /* y.U*Q = b  (if possible) */
    ok = spasm_dense_forward_solve(U, y, SPASM_IDENTITY_PERMUTATION);

    if (ok == SPASM_SUCCESS) {
      /* y.LUQ = b */
      spasm_dense_back_solve(L, y);

      /* x.PLUQ = b */
      spasm_ipvec(PLUQ->p, y, x, n);
    }

    free(y) ;
    spasm_free_LU(PLUQ);
    return ok;
}
