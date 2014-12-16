#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *B, *C;
  int n, m, test, i, j, px, k, nb;
  int *rr, *p, *x, *pinv, *Cp, *Cj;
  spasm_partition *P;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  assert(n == m);

  // generate random row & col permutation
  p = spasm_random_permutation(n);
  pinv = spasm_pinv(p, n);
  B = spasm_permute(A, p, pinv, SPASM_WITH_NUMERICAL_VALUES);

  free(pinv);
  spasm_csr_free(A);

  P = spasm_strongly_connected_components(B);
  p = P->p;
  rr = P->rr;
  nb = P->nr;

  /* verbosity
  printf("p = ");
  for(k = 0; k < n; k++) {
    printf("%d ", p[k] + 1);
  }
  printf("\n -----------------------\n");
  printf("r = ");
  for(k = 0; k < nb + 1; k++) {
    printf("%d ", r[k]);
  }
  printf("\n ");
  */

  /* --- check that p is actually a permutation ---------------- */
  x = spasm_malloc(n * sizeof(int));

  for(i = 0; i < n; i++) {
    x[i] = 0;
  }
  for(i = 0; i < n; i++) {
    if (p[i] < 0 || p[i] >= n) {
      printf("not ok %d - SCC - p is out of range, p[%d] = %d\n", test, i, p[i]);
      exit(0);
    }
    x[ p[i] ]++;
  }

  for(i = 0; i < n; i++) {
    if (x[i] != 1) {
      printf("not ok %d - SCC - p is not bijective\n", test);
      exit(0);
    }
  }

  free(x);

  /* --- verbosity ---------------- */
  printf("#SCC = %d\n", nb);

  for(k = 0; k < nb; k++) {
    printf("# SCC_%d : ", k);
    for(i = rr[k]; i < rr[k + 1]; i++) {
      printf("%d ", p[i]);
    }
    printf("\n");
  }

  /* --- check that decomposition is really block-upper-triangular ---------------- */
  pinv = spasm_pinv(p, n);
  C = spasm_permute(B, p, pinv, SPASM_IGNORE_VALUES);
  Cp = C->p;
  Cj = C->j;

  free(pinv);
  spasm_csr_free(B);

  for(k = 0; k < nb; k++) {
    for(i = rr[k]; i < rr[k + 1]; i++) {
      for(px = Cp[i]; px < Cp[i + 1]; px++) {
	j = Cj[px];
	if (j < rr[k]) {
	  printf("not ok %d - SCC - row %d (in C_%d) has entries on column %d\n", test, i, k, j);
	  exit(0);
	}
      }
    }
  }

  printf("ok %d - SCC\n", test);

  spasm_csr_free(C);
  return 0;
}
