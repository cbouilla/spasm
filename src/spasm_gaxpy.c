#include <assert.h>
#include "spasm.h"

/*
 * (sparse) Matrix x (dense vector)
 */

/* y <--- y + A*x
 */
void spasm_gaxpy(const spasm * A, const spasm_GFp *x, spasm_GFp *y) {
  int j, k, n, prime;
  int *Ap, *Ai, *Ax;

    /* check inputs */
    assert(x != NULL);
    assert(y != NULL);
    assert( spasm_is_csc(A) );

    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    prime = A->prime;

    for (j = 0 ; j < n ; j++) {
      spasm_scatter(Ai, Ax, Ap[j], Ap[j + 1], x[j], y, prime);
    }
}
