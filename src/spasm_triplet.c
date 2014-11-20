#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

/* add an entry to a triplet matrix; enlarge it if necessary */
void spasm_add_entry(spasm *T, int i, int j, spasm_GFp x) {
  assert(spasm_is_triplet(T));
  assert((i >= 0) && (j >= 0));

  if (T->nz == T->nzmax) {
    spasm_sprealloc(T, 2 * T->nzmax);
  }

  if (T->x != NULL) {
    T->x[T->nz] = x;
  }
  T->i[T->nz] = i;
  T->p[T->nz++] = j;
  T->m = spasm_max(T->m, i + 1);
  T->n = spasm_max(T->n, j + 1);
}

/* C = compressed-column form of a triplet matrix T */
spasm * spasm_compress(const spasm *T) {
  int m, n, nz, sum, p, k, *Cp, *Ci, *w, *Ti, *Tj;
    spasm_GFp *Cx, *Tx;
    spasm *C;
    assert(spasm_is_triplet(T));

    m = T->m;
    n = T->n;
    Ti = T->i;
    Tj = T->p;
    Tx = T->x;
    nz = T->nz;
    C = spasm_spalloc(m, n, nz, T->prime, Tx != NULL, false);    /* allocate result */
    w = spasm_calloc(n, sizeof(int));                    /* get workspace */
    Cp = C->p;
    Ci = C->i;
    Cx = C->x;

    /* compute column counts */
    for (k = 0; k < nz; k++) {
        w[Tj[k]]++;
    }
    /* compute column pointers (in both Cp and w) */
    sum = 0;
    for(k = 0; k < n; k++) {
      Cp[k] = sum;
      sum += w[k];
      w[k] = Cp[k];
    }
    Cp[n] = sum;
    /* dispatch entries */
    for (k = 0; k < nz; k++) {
      p = w[Tj[k]]++;            /* A(i,j) is the p-th entry in C */
      Ci[p] = Ti[k];
      if (Cx != NULL) {
	Cx[p] = Tx[k];
      }
    }
    /* success; free w and return C */
    free(w);
    return C;
}
