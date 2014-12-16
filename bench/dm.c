#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *B, *H, *S, *V;
  spasm_partition *DM, *P;
  int n, m, h, s, v, k, largest, ns, i;
  int *rr, *Prr, *cc, *p, *q, *qinv;


  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  printf("A : %d x %d with %d nnz\n", n, m, spasm_nnz(A));

  // compute DM decomposition of A.
  DM = spasm_dulmage_mendelson(A);
  p = DM->p;
  q = DM->q;
  rr = DM->rr;
  cc = DM->cc;

  printf("rr = %d %d %d %d %d\n", rr[0], rr[1], rr[2], rr[3], rr[4]);
  printf("cc = %d %d %d %d %d\n", cc[0], cc[1], cc[2], cc[3], cc[4]);

  h = (cc[2] != 0);
  s = (rr[2] != rr[1]);
  v = (rr[4] != rr[2]);

  qinv = spasm_pinv(q, m);
  B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);
  spasm_csr_free(A);
  free(qinv);

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
    Prr = P->rr;
    k = P->nr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, Prr[i + 1] - Prr[i]);
      ns += (Prr[i + 1] - Prr[i] > 1);
    }
    if (k > 1) {
      printf("H' (%d x %d): %d strongly connected components, %d non-singleton, largest = %.1f %%\n", H->n, H->m, k, ns, 100.0 * largest / H->n);
    }
    spasm_csr_free(H);
    spasm_partition_free(P);
  }

  if (s) {
    S = spasm_submatrix(B, rr[1], rr[2], cc[2], cc[3], SPASM_IGNORE_VALUES);
    printf("S (%d x %d -- %d nnz)\n", S->n, S->m, spasm_nnz(S));
    assert(S->n == S->m);

    P = spasm_strongly_connected_components(H);
    printf("got back\n");
    Prr = P->rr;
    k = P->nr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, Prr[i + 1] - Prr[i]);
      ns += (Prr[i + 1] - Prr[i] > 1);
    }
    if (k > 0) {
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
    Prr = P->rr;
    k = P->nr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, Prr[i + 1] - Prr[i]);
      ns += (Prr[i + 1] - Prr[i] > 1);
    }
    if (k > 1) {
      printf("V' (%d x %d) : %d strongly connected components, %d non-singleton, largest = %.1f %%\n", V->n, V->m, k, ns, 100.0 * largest / V->n);
    }
    spasm_csr_free(V);
    spasm_partition_free(P);
  }

  spasm_partition_free(DM);
  spasm_csr_free(B);
  return 0;
}
