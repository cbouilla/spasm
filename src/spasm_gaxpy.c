#include <assert.h>
#include "spasm.h"

/*
 * (sparse) Matrix x (dense vector)
 */

/* y <--- y + A*x
 */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp *y) {
  int p, j, n, prime, *Ap, *Ai ;
    spasm_GFp *Ax ;

    /* check inputs */
    assert(x != NULL);
    assert(y != NULL);
    assert(spasm_is_csc(A));

    n = A->n ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    prime = A->prime;

    for (j = 0 ; j < n ; j++) {
        for (p = Ap [j] ; p < Ap [j+1] ; p++) {
	  /* TODO : this is axpy-inplace */
	  y [Ai [p]] = (y[Ai[p]] + (Ax[p] * x[j] % prime)) % prime ;
        }
    }
}
