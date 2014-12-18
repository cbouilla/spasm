#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *B, *C, *B_t;
  spasm_partition *DM;
  int n, m, test, i, j, px;
  int *rr, *cc, *p, *q, *x, *y, *Cp, *Cj, *pinv, *qinv, *imatch, *jmatch,
    *Cjmatch, *Cimatch, *Hjmatch, *Sjmatch, *Simatch, *Vimatch;

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

  B = spasm_permute(A, p, q, SPASM_IGNORE_VALUES);
  free(p);
  free(q);
  spasm_csr_free(A);

  B_t = spasm_transpose(B, SPASM_IGNORE_VALUES);
  jmatch = spasm_malloc(n * sizeof(int));
  imatch = spasm_malloc(m * sizeof(int));

   /* --- Maximum matching ------------------------------------------------- */
   if (n < m) {
     spasm_maximum_matching(B, jmatch, imatch);
   } else {
     spasm_maximum_matching(B_t, imatch, jmatch);
   }

  // compute DM decomposition of permuted M.
   DM = spasm_dulmage_mendelson(B, B_t, jmatch, imatch);
   rr = DM->rr;
   cc = DM->cc;
   p = DM->p;
   q = DM->q;

  /* --- check that p and q are actually permutations ---------------- */
  x = spasm_malloc(n * sizeof(int));
  y = spasm_malloc(m * sizeof(int));

  spasm_zero_vector(n, x);
  for(i = 0; i < n; i++) {
    x[ p[i] ]++;
  }
  for(i = 0; i < n; i++) {
    if (x[i] != 1) {
      printf("not ok %d - DM(A) - p is not bijective\n", test);
      exit(0);
    }
  }

  spasm_zero_vector(m, y);
  for(i = 0; i < m; i++) {
    y[ q[i] ]++;
  }
  for(i = 0; i < m; i++) {
    if (y[i] != 1) {
      printf("not ok %d - DM(A) - q is not bijective\n", test);
      exit(0);
    }
  }

  /* --- verbosity ---------------- */
  printf("# sizes : %d x %d | %d x %d | %d x %d\n", rr[1], cc[2], rr[2] - rr[1], cc[3] - cc[2], rr[4] - rr[2], cc[4] - cc[3]);

  /*  for(j = 1; j < 4; j++) {
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
    }*/

  /* --- check that coarse decomposition is really block-upper-triangular ---------------- */
  pinv = spasm_pinv(p, n);
  qinv = spasm_pinv(q, m);
  C = spasm_permute(B, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  Cp = C->p;
  Cj = C->j;
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

  /* --- permute the matching and check the result ------------------------- */
  Cjmatch = spasm_permute_row_matching(n, jmatch, p, qinv);
  Cimatch = spasm_permute_column_matching(m, imatch, pinv, q);

  // check that it is still bijective
  spasm_zero_vector(m, y);
  for(i = 0; i < n; i++) {
    if (Cjmatch[i] != -1) {
      y[ Cjmatch[i] ]++;
    }
  }
  for(j = 0; j < m; j++) {
    if (y[j] > 1) {
      printf("not ok %d - DM - permuted row matching no longer bijective\n", test);
      exit(0);
    }
  }

  spasm_zero_vector(n, x);
  for(j = 0; j < m; j++) {
    if (Cimatch[j] != -1) {
      x[ Cimatch[j] ]++;
    }
  }
  for(i = 0; i < n; i++) {
    if (x[i] > 1) {
      printf("not ok %d - DM - permuted column matching no longer bijective\n", test);
      exit(0);
    }
  }

  // check that the matching is correct (i.e. the right entries exist in the matrix)
  for(i = 0; i < n; i++) {
    // skip unmatched rows
    if (Cjmatch[i] == -1) {
      continue;
    }

    j = 0;
    for(px = Cp[i]; px < Cp[i + 1]; px++) {
      if (Cj[px] == Cjmatch[i]) {
	j = 1;
	break;
      }
    }
    if (j == 0) {
      printf("not ok %d - DM - permuted matching does not follow C\n", test);
      exit(0);
    }
  }

  // check that matching restricted to H is row-perfect.
  Hjmatch = spasm_submatching(Cjmatch, rr[0], rr[1]);
  for(i = rr[0]; i < rr[1]; i++) {
    if (Hjmatch[ i - rr[0] ] == -1) {
      printf("not ok %d - DM - permuted row matching no longer row-perfect on H\n", test);
      exit(0);
    }
  }
  free(Hjmatch);


  // check that matching restricted to S is perfect (both row-perfect and column-perfect).
  Sjmatch = spasm_submatching(Cjmatch, rr[1], rr[2]);
  Simatch = spasm_submatching(Cimatch, cc[2], cc[3]);
  for(i = rr[1]; i < rr[2]; i++) {
    if (Sjmatch[ i - rr[1] ] == -1) {
      printf("not ok %d - DM - permuted row matching no longer row-perfect on S\n", test);
      exit(0);
    }
  }

  for(j = cc[2]; j < cc[3]; j++) {
    if (Simatch[ j - cc[2] ] == -1) {
      printf("not ok %d - DM - permuted column matching no longer column-perfect on S\n", test);
      exit(0);
    }
  }
  free(Sjmatch);
  free(Simatch);

  // check that matching restricted to V is column-perfect.
  Vimatch = spasm_submatching(Cimatch, cc[3], cc[4]);
  for(j = cc[3]; j < cc[4]; j++) {
    if (Vimatch[ j - cc[3] ] == -1) {
      printf("not ok %d - DM - permuted column matching no longer column-perfect on V\n", test);
      exit(0);
    }
  }
  free(Vimatch);


  printf("ok %d - DM(A)\n", test);

  free(Cjmatch);
  free(Cimatch);
  free(x);
  free(y);
  free(qinv);
  spasm_csr_free(C);
  spasm_partition_free(DM);
  return 0;
}
