#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

/* add an entry to a triplet matrix; enlarge it if necessary */
void spasm_add_entry(spasm_triplet *T, int i, int j, spasm_GFp x) {
  int prime;

  assert(T != NULL);
  assert((i >= 0) && (j >= 0));

  prime = T->prime;

  if (T->nz == T->nzmax) {
    spasm_triplet_realloc(T, 2 * T->nzmax);
  }

  if (T->x != NULL) {
    T->x[ T->nz ] = ((x % prime) + prime) % prime;
  }
  T->i[ T->nz ] = i;
  T->j[ T->nz ] = j;
  T->nz += 1;
  T->n = spasm_max(T->n, i + 1);
  T->m = spasm_max(T->m, j + 1);
}

/* C = compressed-row form of a triplet matrix T */
spasm * spasm_compress(const spasm_triplet *T) {
  int m, n, nz, sum, p, k, *Cp, *Cj, *w, *Ti, *Tj;
    spasm_GFp *Cx, *Tx;
    spasm *C;

    m = T->m;
    n = T->n;
    Ti = T->i;
    Tj = T->j;
    Tx = T->x;
    nz = T->nz;

    /* allocate result */
    C = spasm_csr_alloc(n, m, nz, T->prime, Tx != NULL);

    /* get workspace */
    w = spasm_calloc(n, sizeof(int));
    Cp = C->p;
    Cj = C->j;
    Cx = C->x;

    /* compute row counts */
    for (k = 0; k < nz; k++) {
	w[ Ti[k] ]++;
    }

    /* compute row pointers (in both Cp and w) */
    sum = 0;
    for(k = 0; k < n; k++) {
      Cp[k] = sum;
      sum += w[k];
      w[k] = Cp[k];
    }
    Cp[n] = sum;

    /* dispatch entries */
    for (k = 0; k < nz; k++) {
      p = w[ Ti[k] ]++;            /* A(i,j) is the p-th entry in C */
      Cj[p] = Tj[k];
      if (Cx != NULL) {
	Cx[p] = Tx[k];
      }
    }

    /* success; free w and return C */
    free(w);
    return C;
}
