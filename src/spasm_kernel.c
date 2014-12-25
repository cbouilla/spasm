#include <assert.h>
#include "spasm.h"

/* Returns a basis of the (left) kernel of A
 */
spasm * spasm_kernel(const spasm *A) {
    spasm_lu *PLUQ;
    spasm *A_t, *L, *K;
    int i, j, m, r, top, prime, nz;
    int *row_permutation, *xi, *Kp, *Kj, *Kx, *q;
    spasm_GFp *x; //, *y;

    /* check inputs */
    assert(A != NULL);
    A_t = spasm_transpose(A, SPASM_WITH_NUMERICAL_VALUES);

    row_permutation = spasm_cheap_pivots(A_t);
    PLUQ = spasm_PLUQ(A_t, row_permutation, SPASM_DISCARD_L);
    spasm_csr_free(A_t);
    free(row_permutation);

    L = spasm_transpose(PLUQ->U, SPASM_WITH_NUMERICAL_VALUES);

    //    assert(spasm_is_lower_triangular(L));

    m = L->n;
    r = L->m;
    prime = L->prime;
    q = spasm_pinv(PLUQ->qinv, m);

    /* allocate result and workspace */
    K = spasm_csr_alloc(m - r, m, L->nzmax, prime, SPASM_WITH_NUMERICAL_VALUES);
    xi = malloc(3*r * sizeof(int));
    x = malloc(m * sizeof(spasm_GFp));
    //    y = malloc(r * sizeof(spasm_GFp));

    nz = 0;
    Kp = K->p;
    for(i = r; i < m; i++) {
      //      spasm_vector_zero(x, m);
      //      spasm_vector_zero(y, r);
      top = spasm_sparse_backward_solve(L, L, i, xi, x, SPASM_IDENTITY_PERMUTATION);
      //      printf("r = %d, top = %d\n", r, top);

      /* ---- test ----
      spasm_gaxpy(L, x, y);
      spasm_scatter(L->j, L->x, L->p[i], L->p[i + 1], prime - 1, y, prime);
      for(j = 0; j < r; j++) {
	assert( y[j] == 0 );
      }
      */

      /* enlarge K if necessary */
      if (nz + m - top + 1 > K->nzmax) {
	spasm_csr_realloc(K, 2 * K->nzmax + m - top + 1);
      }
      Kp = K->p;
      Kj = K->j;
      Kx = K->x;

      /* finalize previous row of K */
      Kp[i - r] = nz;

      for(j = top; j < r; j++) {
	Kj[nz] = q[ xi[j] ];
	Kx[nz] = x[ xi[j] ];
	nz++;
      }
      Kj[nz] = q[ i ];
      Kx[nz] = prime - 1;
      nz++;
    }

    /* finalize last row */
    Kp[m - r] = nz;

    /* cleanup and return */
    free(xi);
    free(x);
    free(q);
    spasm_free_LU(PLUQ);
    spasm_csr_free(L);
    return K;
}
