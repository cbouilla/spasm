#include <stdio.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *A_t, *B;
  spasm_partition *DM;
  int x, y, n, m, *p, *q, *qinv, *imatch, *jmatch;

  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  x = m;
  y = n;
  if (argc > 1) {
    x = atoi(argv[1]);
  }
  if (argc > 2) {
    y = atoi(argv[2]);
  }

  A_t = spasm_transpose(A, SPASM_IGNORE_VALUES);
  jmatch = spasm_malloc(n * sizeof(int));
  imatch = spasm_malloc(m * sizeof(int));

   /* --- Maximum matching ------------------------------------------------- */
   if (n < m) {
     spasm_maximum_matching(A, jmatch, imatch);
   } else {
     spasm_maximum_matching(A_t, imatch, jmatch);
   }

  // compute DM decomposition of permuted M.
   DM = spasm_dulmage_mendelson(A, A_t, jmatch, imatch);

  p = DM->p;
  q = DM->q;

  qinv = spasm_pinv(q, m);
  B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);
  spasm_csr_free(A);
  free(qinv);

  spasm_save_ppm(stdout, x, y, B, DM);

  spasm_partition_free(DM);
  spasm_csr_free(B);
  return 0;
}
