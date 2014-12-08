#include <assert.h>
#include "spasm.h"

/* reads a sparse matrix from file f in "coordinate text file".
 * (cf. http://math.nist.gov/MatrixMarket/formats.html)
 *
 * return a triplet form matrix.
 *
spasm_triplet * spasm_load_triplet(FILE *f, int prime) {
    int i, j;
    spasm_GFp x ;
    spasm_triplet *T ;
    assert(f != NULL);

    T = spasm_triplet_alloc (0, 0, 1, prime, 1);

    while (fscanf (f, "%d %d %d\n", &i, &j, &x) == 3) {
      spasm_add_entry(T, i, j, x);
    }

    return T;
}
*/

spasm_triplet * spasm_load_sms(FILE *f, int prime) {
    int i, j ;   /* use double for integers to avoid int conflicts */
    spasm_GFp x ;
    spasm_triplet *T;
    char type;
    assert(f != NULL);

    if (fscanf (f, "%d %d %c\n", &i, &j, &type) != 3) {
      fprintf(stderr, "[spasm_load_sms] bad SMS file (header)\n");
      exit(0);
    }

    if (type != 'M') {
      fprintf(stderr, "[spasm_load_sms] only ``Modular'' type supported\n");
      exit(0);
    }

    /* allocate result */
    T = spasm_triplet_alloc (i, j, 1, prime, 1);

    while (fscanf (f, "%d %d %d\n", &i, &j, &x) == 3) {
      if (i == 0 && j == 0 && x == 0) {
	break;
      }
      assert(i != 0);
      assert(j != 0);
      spasm_add_entry(T, i - 1, j - 1, x);
    }

    return T;
}


void spasm_save_csr(FILE *f, const spasm *A) {
  int i, n, p;
  int *Aj, *Ap;
  spasm_GFp *Ax;

    assert(f != NULL);
    assert(A != NULL);

    Aj = A->j;
    Ap = A->p;
    Ax = A->x;
    n  = A->n;

    /* compressed row form */
    for(i = 0; i < n; i++) {
      for(p = Ap[i]; p < Ap[i + 1]; p++) {
	fprintf(f, "%d %d %d\n", i, Aj[p], (Ax != NULL) ? Ax[p] : 1);
      }
    }
}

void spasm_save_sms(FILE *f, const spasm_triplet *A) {
  int i, nz, n, m;
  int *Ai, *Aj;
  spasm_GFp *Ax;

    assert(f != NULL);
    assert(A != NULL);
    Ai = A->i;
    Aj = A->j;
    Ax = A->x;
    nz = A->nz;
    n  = A->n;
    m  = A->m;

    fprintf(f, "%d %d M\n", n, m);

    for(i = 0; i < nz; i++) {
      fprintf(f, "%d %d %d\n", Ai[i] + 1, Aj[i] + 1, (Ax != NULL) ? Ax[i] : 1);
    }

    fprintf(f, "0 0 0\n");
}
