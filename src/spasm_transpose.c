#include <assert.h>
#include "spasm.h"

spasm *spasm_transpose(const spasm *C, int keep_values) {
  int i, j, m, n, sum, px, py, *Cp, *Cj, *w, *Tp, *Tj;
  spasm_GFp *Cx, *Tx;
  spasm *T;

  m = C->m;
  n = C->n;
  Cp = C->p;
  Cj = C->j;
  Cx = C->x;

  /* allocate result */
  T = spasm_csr_alloc(m, n, spasm_nnz(C), C->prime, keep_values && (Cx != NULL));
  Tp = T->p;
  Tj = T->j;
  Tx = T->x;

  /* get workspace */
  w = spasm_calloc(m, sizeof(int));

  /* compute column counts */
  for (i = 0; i < n; i++) {
    for(px = Cp[i]; px < Cp[i + 1]; px++) {
      j = Cj[px];
      w[j]++;
    }
  }

  /* compute column pointers (in both Cp and w) */
  sum = 0;
  for(j = 0; j < m; j++) {
    Tp[j] = sum;
    sum += w[j];
    w[j] = Tp[j];
  }
  Tp[m] = sum;

  /* dispatch entries */
  for (i = 0; i < n; i++) {
    for(px = Cp[i]; px < Cp[i + 1]; px++) {
      j = Cj[px];
      py = w[j];
      Tj[py] = i;
      if (Tx != NULL) {
	Tx[py] = Cx[px];
      }
      w[j]++;
    }
  }

  free(w);
  return T;
}
