#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "spasm.h"

/** permute a matrix so that all empty rows are removed. Also generate
    the permutation that reproduces the original.
 */

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A;
  int n, m, i, j, t, u, v;
  int *p, *q, *Ap, *Aj;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  Ap = A->p;
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
  for (i = 0; i < n; i++) {
    for(t = Ap[i]; t < Ap[i + 1]; t++) {
      j = Aj[t];
      p[i] = 1;
      q[j] = 1;
    }
  }

  /* count non-empty rows/columns */
  u = 0;
  for(i = 0; i<n; i++) {
    u += p[i];
  }
  v = 0;
  for(j = 0; j<m; j++) {
    v += q[j];
  }
  
  printf("matrix has advertized dimension %d x %d but is in fact %d x %d\n", n, m, u, v);
  
  spasm_csr_free(A);
  free(p);
  free(q);
  return 0;
}
