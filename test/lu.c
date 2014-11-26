/*
 * sage utility script to generate random sparse matrix
 *
sage: F = GF(257)
sage: n = 100
sage: out = open("mat.txt", "w")
sage: M = random_matrix(F, n, n, density=0.25, sparse=True)
sage: for (i,j) in M.nonzero_positions():
....:     out.write("{0} {1} {2}\n".format(i, j, M[i,j]))
....:
....: out.close()
sage: V = F^n
sage: x = V([i + 1 for i in range(n)])
sage: x*M
*/
#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *C;
  spasm_lu *LU;
  int i, n;

  T = spasm_load_triplet(stdin, 257);
  C = spasm_compress(T);
  spasm_triplet_free(T);

  LU = spasm_LU(C);

  spasm_save_csr(stdout, LU->L);
  fprintf(stdout, "--------------------------------------\n");
  spasm_save_csr(stdout, LU->U);

  spasm_csr_free(C);
  return 0;
}
