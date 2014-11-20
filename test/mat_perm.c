#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  int n, m, i, test;
  spasm *T, *A, *C;
  int *pinv, *q;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_ctf(stdin, 257);
  A = spasm_compress(T);
  spasm_spfree(T);
  n = A->n;
  m = A->m;

  pinv = malloc(m * sizeof(int));
  for(i = 0; i < m; i++) {
    pinv[i] = m - i - 1;
  }
  q = malloc(n * sizeof(int));
  for(i = 0; i < n; i++) {
    q[i] = n - i - 1;
  }

  switch(test) {
  // test 1 : only reverse line order
  case 1:
    C = spasm_permute(A, pinv, SPASM_IDENTITY_PERMUTATION, 1);
    break;

    // test 1 : reverse both line and column order (this should transpose)
  case 2:
    C = spasm_permute(A, pinv, q, 1);
    break;
  };

  spasm_spfree(A);
  spasm_save_ctf(stdout, C);
  spasm_spfree(C);
}
