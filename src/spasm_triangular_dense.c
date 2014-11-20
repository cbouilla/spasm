#include <assert.h>
#include "spasm.h"

/*
 * Solving triangular systems, dense RHS
 */

/* solve L.x=b where x and b are dense.  x=b on input, solution on output.
 *
 * Assumption : the diagonal entry is always present, is always != 0, and is the first one of each column.
 * More precisely, L[j,j] is Lx[ Lp[j] ]
 */
void spasm_lsolve(const spasm * L, spasm_GFp * x) {
    int p, j, n, *Lp, *Li, prime;
    spasm_GFp *Lx;

    /* check inputs */
    assert(x != NULL);
    assert(spasm_is_csc(L));

    n = L->n;
    Lp = L->p;
    Li = L->i;
    Lx = L->x;
    prime = L->prime;

    for (j = 0; j < n; j++) {
      /* check diagonal entry */
      assert( Li[ Lp[j] ] == j );
      const spasm_GFp diagonal_entry = Lx[Lp[j]];

      x[j] = (x[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
      for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
	/* TODO : this is axpy-inplace */
	x[Li[p]] = (x[Li[p]] + prime - ((Lx[p] * x[j]) % prime)) % prime;
      }
    }
}

/* solve U.x=b where x and b are dense.  x=b on input, solution on output.
 *
 * Assumption : the diagonal entry is always present, is always != 0, and is the last one of each column.
 * More precisely, U[j,j] is Ux[ Up[j+1] - 1 ]
 */
void spasm_usolve (const spasm *U, spasm_GFp *x) {
  int p, j, n, *Up, *Ui, prime ;
  spasm_GFp *Ux ;

    /* check inputs */
    assert(x != NULL);
    assert(spasm_is_csc(U));

    n = U->n ;
    Up = U->p ;
    Ui = U->i ;
    Ux = U->x ;
    prime = U->prime;
    
    for (j = n-1 ; j >= 0 ; j--) {
      /* check diagonal entry */
      printf("Column %d, diagonal entry row = %d\n", j, Ui[ Up[j + 1] - 1 ]);
      assert( Ui[ Up[j + 1] - 1 ] == j );
      const spasm_GFp diagonal_entry = Ux[ Up[ j + 1 ] - 1 ];

        x[j] = (x[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
        for (p = Up [j] ; p < Up [j+1]-1 ; p++) {
            x[ Ui[ p ] ] = (x[ Ui[p] ] + prime - ((Ux[p] * x[j]) % prime)) % prime;
        }
    }
}
