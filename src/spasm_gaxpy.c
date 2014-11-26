#include <assert.h>
#include "spasm.h"

/*
 * (sparse) Matrix x (dense vector)
 */

/* y <--- y + x*A
 */
void spasm_gaxpy(const spasm * A, const spasm_GFp *x, spasm_GFp *y) {
  int i, n, prime;
  int *Ap, *Aj, *Ax;

    /* check inputs */
    assert(x != NULL);
    assert(y != NULL);
    assert(A != NULL);

    n = A->n;
    Ap = A->p;
    Aj = A->j;
    Ax = A->x;
    prime = A->prime;

    for (i = 0 ; i < n ; i++) {
      spasm_scatter(Aj, Ax, Ap[i], Ap[i + 1], x[i], y, prime);
    }
}
