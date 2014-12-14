#include <assert.h>
#include "spasm.h"

/* returns the number of non-trivial strongly connected component of A
   (which must be square), seen as an directed graph. */
int spasm_strongly_connected_components(const spasm *A) {
  int n, m, i, k, size, n_scc, top;
  int *pstack, *mark, *xi, *p;
  spasm *A_t;

  n = A->n;
  m = A->m;
  assert(n == m);

  pstack = spasm_malloc(n * sizeof(int));
  mark = spasm_malloc(n * sizeof(int));
  xi = spasm_malloc(n * sizeof(int));
  p = spasm_malloc(n * sizeof(int));

  // first pass
  for (i = 0 ; i < n ; i++) {
    mark[i] = 0;
  }

  for (i = 0 ; i < n ; i++) {
    if (mark[i]) {
      continue;
    }
    top = spasm_dfs(i, A, top, xi, pstack, mark, SPASM_IDENTITY_PERMUTATION);
  }
  assert(top == 0);

  // second pass : in the order given by xi, in A_t
  for (i = 0 ; i < n ; i++) {
    mark[i] = 0;
  }

  A_t = spasm_transpose(A, SPASM_IGNORE_VALUES);
  top = n;

  //  nb = n;
  n_scc = 0;
  for(k = 0; k < n; k++) {
    i = xi[k];

    if (mark[i]) {
      continue;
    }
    // r[nb] = top;
    // nb--;
    spasm_dfs(i, A_t, top, p, pstack, mark, SPASM_IDENTITY_PERMUTATION);
    n_scc++;
  }

  free(xi);
  free(pstack);
  free(mark);
  spasm_csr_free(A_t);
  return n_scc;
}
