/*
 * sage utility script to generate random sparse matrix
 *
sage: F = GF(7)
sage: n = 128
sage: M = matrix(F, n, n, sparse=True)
sage: while M.is_singular():
sage:     M = random_matrix(F, n, n, density=0.33, sparse=True)
sage: out = open("mat.txt", "w")
sage: for (i,j) in M.nonzero_positions():
....:     out.write("{0} {1} {2}\n".format(i, j, M[i,j]))
....:
sage: out.close()
*/
#include <stdio.h>
#include <assert.h>
#include "spasm.h"

/*
[2 0 0 4]
[2 0 6 0]
[5 0 6 0]
[0 6 0 0]
*/

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A;
  spasm_lu *LU;
  spasm_GFp *x, *y, *u, *v;
  int n, m, test, i, j, fail;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_triplet(stdin, 7);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(m * sizeof(spasm_GFp));
  u = malloc(n * sizeof(spasm_GFp));
  v = malloc(m * sizeof(spasm_GFp));

  LU = spasm_LU(A);

  for(i = 0; i < n; i++) {
    printf("checking i = %d\n", i);

    for(j = 0; j < n; j++) {
      x[j] = 0;
      u[j] = 0;
    }
    for(j = 0; j < m; j++) {
      y[j] = 0;
      v[j] = 0;
    }
    x[i] = 1;

    // ce test est-il correct ?
    spasm_gaxpy(A, x, y); // y <- x*A
    spasm_gaxpy(LU->L, x, u); // u <- x*L
    spasm_gaxpy(LU->U, u, v); // v <- (x*L)*U

    fail = 0;
    for(j = 0; j < m; j++) {
      if (y[j] != v[j]) {
	printf("not ok %d - L*U == A (col %d)\n", test, j);
	exit(1);
      }
    }
  }

  printf("ok %d - L*U == A\n", test);

  spasm_csr_free(A);
  spasm_free_LU(LU);
  free(x);
  free(y);
  free(u);
  free(v);
  
  return 0;
}
