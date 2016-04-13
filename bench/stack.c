#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "spasm.h"

/** permute a matrix so that all empty rows are removed. Also generate
    the permutation that reproduces the original.
 */

int main() {
  spasm_triplet *A;
  int n, m, nz, i, j, k, u, v;
  int *p, *q, *Ai, *Aj;

  A = spasm_load_sms(stdin, 42013);
  
  n = A->n;
  m = A->m;
  nz = A->nz;
  Ai = A->i;
  Aj = A->j;
  
  /* allocate result */
  p = spasm_malloc(n * sizeof(int));
  q = spasm_malloc(m * sizeof(int));
  for(i = 0; i<n; i++) {
    p[i] = 0;
  }
  for(j = 0; j<m; j++) {
    q[j] = 0;
  }

  /* mark non-empty rows/columns */
  for (k = 0; k < nz; k++) {
      p[ Ai[k] ] = 1;
      q[ Aj[k] ] = 1;
  }

  /* sum-prefix: p[i] = sum(p[k], k in range(i)) */
  u = 0;
  for(i = 0; i<n; i++) {
    if (p[i] > 0) {
      p[i] = u++;
    }
  }

  v = 0;
  for(j = 0; j<m; j++) {
    if (q[j] > 0) {
      q[j] = v++;
    }
  }
  
  fprintf(stderr, "matrix has advertized dimension %d x %d but is in fact %d x %d\n", n, m, u, v);
  
  /* modify matrix */
  A->n = u;
  A->m = v;
  for (k = 0; k < nz; k++) {
      Ai[k] = p[ Ai[k] ];
      Aj[k] = q[ Aj[k] ];
  }


  spasm_save_triplet(stdout, A);
  spasm_triplet_free(A);
  free(p);
  free(q);
  return 0;
}
