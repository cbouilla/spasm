#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *B;
  int i, j, k, n, m, r;
  int *p, *qinv, *imatch, *jmatch;
  double start_time;
  int transposed = 0;

  T = spasm_load_sms(stdin, 42013);
  n = T->n;
  m = T->m;
  fprintf(stderr, "A : %d x %d with %d nnz (density = %.3f %%)\n", T->n, T->m, T->nz, 100.0 * T->nz / (1.0 * T->n * T->m));
  
  if (T->n > T->m) {
    fprintf(stderr, "[matching - info] transposing matrix to speed-up maximum matching\n");
    spasm_triplet_transpose(T);
    transposed = 1;
  }
  A = spasm_compress(T);

  imatch = spasm_malloc(m * sizeof(int));
  jmatch = spasm_malloc(n * sizeof(int));

  fprintf(stderr, "[matching] maximum matching...");
  fflush(stderr);
  start_time = spasm_wtime();
  
  r = spasm_maximum_matching(A, jmatch, imatch);

  if (transposed) {
    spasm_csr_free(A);
    spasm_triplet_transpose(T);
    A = spasm_compress(T);
    
    p = imatch;
    imatch = jmatch;
    jmatch = p;
  }

  fprintf(stderr, " %.1f s\n", spasm_wtime() - start_time);
  if (r < spasm_min(n, m)) {
    fprintf(stderr, "[matching] matrix is structurally rank-deficient. Structural rank: %d\n", r);
  } else if (n == m) {
    fprintf(stderr, "[matching] matrix is square and structurally full-rank\n");
  } else if (n < m) {
    fprintf(stderr, "[matching] matrix has full structural row-rank\n");  
  } else {
    fprintf(stderr, "[matching] matrix has full structural column-rank\n");  
  }
  
  /* build the permutations */
  p = spasm_malloc(n * sizeof(int));
  qinv = spasm_malloc(m * sizeof(int));
  k = 0;

  /* put matched rows and columns first, in-order */
  for(i = 0; i < n; i++) {
    j = jmatch[i];
    if (j == -1) { /* row i is unmatched */
      continue;
    }
    p[k] = i;
    qinv[j] = k;
    k++;
  }
  assert (k == r);

  /* add unmatched rows */
  for(i = 0; i < n; i++) {
    if (jmatch[i] == -1) {
      p[k] = i;
      k++;
    }
  }
  assert(k == n);

  /* add unmatched columnss */
  k = r;
  for(j = 0; j < m; j++) {
    if (imatch[j] == -1) {
      qinv[j] = k;
      k++;
    }
  }
  assert(k == m);

  /* print permuted matrix */
  B = spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  spasm_save_csr(stdout, B);

  spasm_csr_free(B);
  spasm_csr_free(A);
  free(p);
  free(qinv);
  free(imatch);
  free(jmatch);

  return 0;
}