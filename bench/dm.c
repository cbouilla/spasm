#include <stdio.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *B, *H, *S, *V;
  spasm_partition *DM, *P;
  int n, m, h, s, v, k, largest, ns, i;
  int *rr, *cc, *p, *q, *qinv, *rmark, *cmark, *r;


  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  // compute DM decomposition of A.
  DM = spasm_dulmage_mendelson(A);
  p = DM->p;
  q = DM->q;
  rr = DM->rr;
  cc = DM->cc;

  // --- verbosity ----------------
  k = 0;
  h = (cc[2] != 0);
  s = (rr[2] != rr[1]);
  v = (rr[4] != rr[2]);

  qinv = spasm_pinv(q, m);
  B = spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  spasm_csr_free(A);
  free(qinv);

  rmark = spasm_malloc(n * sizeof(int));
  cmark = spasm_malloc(m * sizeof(int));
  r = spasm_malloc((n+1) * sizeof(int));

  if (h) {
    H = spasm_submatrix(B, rr[0], rr[1], cc[0], cc[2], SPASM_IGNORE_VALUES);
    P = spasm_connected_components(H);
    if (P->nr >= 0) {
      printf("H (%d x %d -- %d nnz): %d connected components\n", H->n, H->m, spasm_nnz(H), P->nr);
    }
    spasm_csr_free(H);
    spasm_partition_free(P);

    H = spasm_submatrix(B, rr[0], rr[1], cc[1], cc[2], SPASM_IGNORE_VALUES);
    P = spasm_strongly_connected_components(H);
    rr = P->rr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, rr[i + 1] - rr[i]);
      ns += (rr[i + 1] - rr[i] > 1);
    }
    if (k > 1) {
      printf("H' (%d x %d): %d strongly connected components, %d non-singleton, largest = %.1f %%\n", H->n, H->m, k, ns, 100.0 * largest / H->n);
    }
    spasm_csr_free(H);
    spasm_partition_free(P);
  }

  if (s) {
    S = spasm_submatrix(B, rr[1], rr[2], cc[2], cc[3], SPASM_IGNORE_VALUES);

    P = spasm_strongly_connected_components(H);
    rr = P->rr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, rr[i + 1] - rr[i]);
      ns += (rr[i + 1] - rr[i] > 1);
    }
    if (k > 1) {
      printf("S (%d x %d) : %d strongly connected components, %d non-singleton, largest = %.1f %%\n", S->n, S->m, k, ns, 100.0 * largest / S->n);
    }
    spasm_csr_free(S);
    spasm_partition_free(P);
  }

  if (v) {
    V = spasm_submatrix(B, rr[2], rr[4], cc[3], cc[4], SPASM_IGNORE_VALUES);
    P = spasm_connected_components(V);
    if (P->nr > 1) {
      printf("V (%d x %d) : %d connected components\n", V->n, V->m, P->nr);
    }
    spasm_csr_free(V);
    spasm_partition_free(P);

    V = spasm_submatrix(B, rr[2], rr[3], cc[3], cc[4], SPASM_IGNORE_VALUES);
    P = spasm_strongly_connected_components(V);
    rr = P->rr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, rr[i + 1] - rr[i]);
      ns += (rr[i + 1] - rr[i] > 1);
    }
    if (k > 1) {
      printf("V' (%d x %d) : %d strongly connected components, %d non-singleton, largest = %.1f %%\n", V->n, V->m, k, ns, 100.0 * largest / V->n);
    }
    spasm_csr_free(V);
    spasm_partition_free(P);
  }

  free(rmark);
  free(cmark);
  spasm_partition_free(DM);
  spasm_csr_free(B);
  return 0;
}
