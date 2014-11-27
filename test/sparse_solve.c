/*
 * sage utility script to generate random sparse lower-triangular matrix
 *
sage: F = GF(257)
sage: n = 100
sage: M = random_matrix(F, n, n, density=0.25, sparse=True)
sage: for i in range(M.nrows()):
....:     M[i,i] = 1
....:     for j in range(i+1, M.ncols()):
....:         M[i,j] = 0
....:
sage: out = open("mat.txt", "w")
sage: for (i,j) in M.nonzero_positions():
....:     out.write("{0} {1} {2}\n".format(i, j, M[i,j]))
....:
sage: out.close()


-- easily adapted to yield an upper-triangular sparse matrix --

sage: for i in range(M.nrows()):
....:     M[i,i] = 1
....:     for j in range(i):
....:         M[i,j] = 0

*/
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *U, *B;
  int fail, i, n, m, test, top, *xi;
  spasm_GFp *x, *y;

  assert(argc == 2);
  test = atoi(argv[1]);

  // load matrix
  T = spasm_load_triplet(stdin, 257);
  U = spasm_compress(T);
  spasm_triplet_free(T);
  n = U->n;
  m = U->m;

  // load RHS
  T = spasm_triplet_alloc(1, m, 10, 257, true);
  spasm_add_entry(T, 0, 0, 1);
  spasm_add_entry(T, 2, 0, 2);
  spasm_add_entry(T, 4, 0, 3);
  B = spasm_compress(T);
  spasm_triplet_free(T);

  xi = malloc(2*n * sizeof(int));
  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(n * sizeof(spasm_GFp));
  for(i = 0; i < n; i++) {
    x[i] = 0;
    y[i] = 0;
  }

  switch(test) {
  case 3: // testing sparse U-solver
    top = spasm_sparse_forwardsolve(U, B, 0, xi, x, SPASM_IDENTITY_PERMUTATION);
    break;
  }

  spasm_gaxpy(U, x, y);
  spasm_scatter(B->j, B->x, B->p[0], B->p[1], B->prime - 1, y, B->prime);

  fail = 0;
  for(i = 0; i < n; i++) {
    fail |= (y[i] != 0);
  }

  if (fail) {
    printf("not ok %d - sparse triangular U-solve\n", test);
  } else {
    printf("ok %d - sparse triangular U-solve\n", test);
  }

  spasm_csr_free(U);
  spasm_csr_free(B);
  return 0;
}
