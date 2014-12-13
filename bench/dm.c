#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A;
  spasm_dm *DM;
  int n, m, k;
  int *rr, *cc;


  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  // compute DM decomposition of A.
  DM = spasm_dulmage_mendelson(A);
  rr = DM->rr;
  cc = DM->cc;

  // --- verbosity ----------------
  k = 0;
  k += (rr[1] != 0) && (cc[2] != 0);
  k += (rr[2] - rr[1] != 0) && (cc[3] - cc[2] != 0);
  k += (rr[4] - rr[2] != 0) && (cc[4] - cc[3] != 0);
  if (k < 2) {
    printf("Dulmage-Mendelsohn decomposition is trivial\n");
  } else {
    printf("Dulmage-Mendelsohn decomposition : %d x %d | %d x %d | %d x %d\n", rr[1], cc[2], rr[2] - rr[1], cc[3] - cc[2], rr[4] - rr[2], cc[4] - cc[3]);
  }
  spasm_dm_free(DM);

  spasm_csr_free(A);
  return 0;
}
