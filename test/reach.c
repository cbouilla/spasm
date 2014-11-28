#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *U;
  int i, n, m, root, *xi, *pstack, *pinv;

  assert(argc > 1);
  root = atoi(argv[1]);

  T = spasm_load_triplet(stdin, 257);
  U = spasm_compress(T);
  spasm_triplet_free(T);

  n = U->n;
  m = U->m;
  assert( n <= m );

  pstack = malloc(n * sizeof(int));
  xi = malloc(n * sizeof(int));

  pinv = NULL;
  if (n < m) { /* upper-trapezoidal */
    pinv = malloc(m * sizeof(int));
    for(i = 0; i < n; i++) {
      pinv[i] = i;
    }
    for(i = n; i < m; i++) {
      pinv[i] = -1;
    }
  }

  i = spasm_dfs(root, U, n, xi, pstack, pinv);
  for( ; i < n; i++) {
    printf("%d\n", xi[i]);
  }

  spasm_csr_free(U);
  return 0;
}
