#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *B, *C;
  spasm_dm *DM;
  int n, m, test, i, j, px;
  int *rr, *cc, *p, *q, *x, *y, *Cp, *Cj, *qinv;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  // generate random row & col permutation
  p = spasm_random_permutation(n);
  q = spasm_random_permutation(m);

  B = spasm_permute(A, p, q, 1);
  free(p);
  free(q);
  spasm_csr_free(A);

  // compute DM decomposition of permuted M.
  DM = spasm_dulmage_mendelson(B);
  rr = DM->rr;
  cc = DM->cc;
  p = DM->p;
  q = DM->q;

  /* --- check that p and q are actually permutations ---------------- */
  x = spasm_malloc(n * sizeof(int));
  y = spasm_malloc(m * sizeof(int));

  for(i = 0; i < n; i++) {
    x[i] = 0;
  }
  for(i = 0; i < n; i++) {
    x[ p[i] ]++;
  }
  for(i = 0; i < n; i++) {
    if (x[i] != 1) {
      printf("not ok %d - DM(A) - p is not bijective\n", test);
      exit(0);
    }
  }

  for(i = 0; i < m; i++) {
    y[i] = 0;
  }
  for(i = 0; i < m; i++) {
    y[ q[i] ]++;
  }
  for(i = 0; i < m; i++) {
    if (y[i] != 1) {
      printf("not ok %d - DM(A) - q is not bijective\n", test);
      exit(0);
    }
  }

  free(x);
  free(y);

  /* --- verbosity ---------------- */
  for(j = 1; j < 4; j++) {
    printf("# R_%d : ", j);
    for(i = rr[j - 1]; i < rr[j]; i++) {
      printf("%d ", p[i] + 1);
    }
    printf("\n");
  }
  printf("# R_0 : ");
  for(i = rr[3]; i < rr[4]; i++) {
    printf("%d ", p[i] + 1);
  }
  printf("\n");

  for(j = 0; j < 4; j++) {
    printf("# C_%d : ", j);
    for(i = cc[j]; i < cc[j + 1]; i++) {
      printf("%d ", q[i] + 1);
    }
    printf("\n");
  }

  /* --- check that coarse decomposition is really block-upper-triangular ---------------- */
  qinv = spasm_pinv(q, m);
  C = spasm_permute(B, p, qinv, 1);
  Cp = C->p;
  Cj = C->j;

  free(qinv);
  spasm_csr_free(B);

  for(i = rr[1]; i < rr[2]; i++) {
    for(px = Cp[i]; px < Cp[i + 1]; px++) {
      j = Cj[px];
      if (j < cc[2]) {
	printf("not ok %d - DM(A) - row %d (in R_2) has entries in C_0 or C_1\n", test, i);
	exit(0);
      }
    }
  }

  for(i = rr[2]; i < rr[4]; i++) {
    for(px = Cp[i]; px < Cp[i + 1]; px++) {
      j = Cj[px];
      if (j < cc[3]) {
	printf("not ok %d - DM(A) - row %d (in R_3 or R_0) has entries in C_0, C_1 or C_2\n", test, i);
	exit(0);
      }
    }
  }

  printf("ok %d - DM(A)\n", test);

  spasm_csr_free(C);
  spasm_dm_free(DM);
  return 0;
}
