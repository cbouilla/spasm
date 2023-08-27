#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"


int main(int argc, char **argv) {
  int n, m, i, j;
  spasm_triplet *T;
  spasm *A, *B, *C;
  spasm_GFp *x, *y, *z;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  B = spasm_transpose(A, SPASM_WITH_NUMERICAL_VALUES);
  C = spasm_transpose(B, SPASM_WITH_NUMERICAL_VALUES);

  if (spasm_is_lower_triangular(A)) {
    printf("# A is lower-triangular.\n");
    printf("Check that A.T is upper-triangular\n");
    if (!spasm_is_upper_triangular(B)) {
      printf("not ok - A.T not upper-triangular \n");
      exit(1);
    }
    printf("Check that A.T.T is lower-triangular\n");
    if (!spasm_is_lower_triangular(C)) {
      printf("not ok - A.T.T not lower-triangular \n");
      exit(1);
    }
  }

  if (spasm_is_upper_triangular(A)) {
    printf("# A is upper-triangular.\n");
    printf("Check that A.T is lower-triangular\n");
    if (!spasm_is_lower_triangular(B)) {
      printf("not ok - A.T not lower-triangular \n");
      exit(1);
    }
    printf("Check that A.T.T is upper-triangular\n");
    if (!spasm_is_upper_triangular(C)) {
      printf("not ok - A.T.T not upper-triangular \n");
      exit(1);
    }
  }

  printf("# check that A.T.T == A\n");
  x = spasm_malloc(n * sizeof(spasm_GFp));
  y = spasm_malloc(m * sizeof(spasm_GFp));
  z = spasm_malloc(m * sizeof(spasm_GFp));

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      x[j] = 0;
    }
    x[i] = 1;
    for(j = 0; j < m; j++) {
      y[j] = 0;
      z[j] = 0;
    }

    spasm_gaxpy(A, x, y);
    spasm_gaxpy(C, x, z);

    for(j = 0; j < m; j++) {
      if (y[j] != z[j]) {
	printf("not ok - (A.T).T != A \n");
	exit(0);
      }
    }
  }
  printf("ok - transpose \n");

  free(x);
  free(y);
  free(z);
  spasm_csr_free(A);
  spasm_csr_free(B);
  spasm_csr_free(C);
  exit(EXIT_SUCCESS);
}