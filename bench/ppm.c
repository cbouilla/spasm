#include <stdio.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *B;
  spasm_partition *DM;
  int x, y, n, m, *p, *q, *qinv;

  T = spasm_load_sms(stdin, -1);
  B = spasm_compress(T);
  spasm_triplet_free(T);
  A = spasm_transpose(B, SPASM_IGNORE_VALUES);
  spasm_csr_free(B);

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

  // compute DM decomposition of A.
  DM = spasm_dulmage_mendelson(A);
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
