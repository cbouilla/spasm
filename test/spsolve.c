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
  spasm *T, *G, *B;
  int fail, i, n, m, test, top, *xi;
  spasm_GFp *x, *y;
  FILE *RHS_file;
  const char * const test_name[]  = {"sparse RHS L-solver", "sparse RHS U-solver"};

  assert(argc == 2);
  test = atoi(argv[1]);

  // load matrix
  T = spasm_load_ctf(stdin, 257);
  G = spasm_compress(T);
  spasm_spfree(T);
  n = G->n;
  m = G->m;

  // load RHS
  T = spasm_spalloc(m, 1, 10, 257, true, true);
  spasm_add_entry(T, 0, 0, 1);
  spasm_add_entry(T, 2, 0, 2);
  spasm_add_entry(T, 4, 0, 3);
  B = spasm_compress(T);
  spasm_spfree(T);

  xi = malloc(2*n * sizeof(int));
  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(n * sizeof(spasm_GFp));
  for(i = 0; i < n; i++) {
    x[i] = 0;
    y[i] = 0;
  }

  switch(test) {
  case 3: // testing sparse L-solver
    top = spasm_sparse_triangular_solve(G, B, 0, xi, x, SPASM_IDENTITY_PERMUTATION, 1);
    break;
  case 4: // testing dense U-solver
    top = spasm_sparse_triangular_solve(G, B, 0, xi, x, SPASM_IDENTITY_PERMUTATION, 0);
    break;
  }

  spasm_gaxpy(G, x, y);
  spasm_scatter(B->i, B->x, B->p[0], B->p[1], B->prime - 1, y, B->prime);

  fail = 0;
  for(i = 0; i < n; i++) {
    fail |= (y[i] != 0);
  }

  if (fail) {
      printf("not ok %d - %s\n", test, test_name[test-3]);
  } else {
      printf("ok %d - %s\n", test, test_name[test-3]);
  }

  return 0;
}
