#include <assert.h>
#include "spasm.h"

/* reads a sparse matrix from file f in "coordinate text file".
 * (cf. http://math.nist.gov/MatrixMarket/formats.html)
 *
 * return a triplet form matrix.
 */
spasm_triplet * spasm_load_triplet(FILE *f, int prime) {
    int i, j ;   /* use double for integers to avoid int conflicts */
    spasm_GFp x ;
    spasm_triplet *T ;
    assert(f != NULL);

    /* allocate result */
    T = spasm_triplet_alloc (0, 0, 1, prime, 1);

    while (fscanf (f, "%d %d %d\n", &i, &j, &x) == 3) {
      spasm_add_entry(T, i, j, x);
    }

    return T;
}

void spasm_save_csr(FILE *f, const spasm *A) {
  int i, n, p;
  int *Ai, *Ap;
  spasm_GFp *Ax;

    assert(f != NULL);
    assert(A != NULL);

    Ai = A->i;
    Ap = A->p;
    Ax = A->x;
    n  = A->n;

    /* compressed column form */
    for(i = 0; i < n; i++) {
      for(p = Ap[i]; p < Ap[i + 1]; p++) {
	fprintf(f, "%d %d %d\n", Ai[p], i, (Ax != NULL) ? Ax[p] : 1);
      }
    }
}

void spasm_save_triplet(FILE *f, const spasm_triplet *A) {
  int i, nz, n;
  int *Ai, *Aj;
  spasm_GFp *Ax;

    assert(f != NULL);
    assert(A != NULL);
    Ai = A->i;
    Aj = A->j;
    Ax = A->x;
    nz = A->nz;
    n  = A->n;

    /* triplet form */
    for(i = 0; i < nz; i++) {
      fprintf(f, "%d %d %d\n", Ai[i], Aj[i], (Ax != NULL) ? Ax[i] : 1);
    }
}
