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
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *G;
  int fail, i, n, test;
  spasm_GFp *x, *y;
  const char * const test_name[]  = {"dense RHS L-solver", "dense RHS U-solver"};

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_triplet(stdin, 257);
  G = spasm_compress(T);
  spasm_triplet_free(T);

  n = G->n;
  assert(n < G->prime);
  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(n * sizeof(spasm_GFp));
  for(i = 0; i < n; i++) {
    x[i] = i + 1;
    y[i] = 0;
  }

  switch(test) {
  case 1: // testing dense L-solver
    spasm_dense_backsolve(G, x);
    break;
  case 2: // testing dense U-solver
    spasm_dense_forwardsolve(G, x);
    break;
  }
// the solution must be in x

  spasm_gaxpy(G, x, y);
  fail = 0;
  for(i = 0; i < n; i++) {
    fail |= (y[i] != i + 1);
  }

  if (fail) {
      printf("not ok %d - %s\n", test, test_name[test-1]);
  } else {
      printf("ok %d - %s\n", test, test_name[test-1]);
  }

  spasm_csr_free(G);
  return 0;
}
