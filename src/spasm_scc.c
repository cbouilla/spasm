/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_scc.c  */
#include <assert.h>
#include "spasm.h"

/*
 * returns the non-trivial strongly connected component of A (which must be
 * square), seen as an directed graph.
 */
/* r must have size n+1 */
spasm_partition *spasm_strongly_connected_components(const spasm * A) {
  int n, m, i, k, n_scc, top;
  int *pstack, *mark, *xi, *rcopy, *p, *q, *rr, *cc;
  spasm *A_t;
  spasm_partition *P;

  n = A->n;
  m = A->m;
  assert(n == m);

  P = spasm_partition_alloc(n, n, n, n);
  /* p and q are identical, but we fill them both anyway */
  p = P->p;
  q = P->q;
  rr = P->rr;
  cc = P->cc;

  pstack = spasm_malloc(n * sizeof(int));
  mark = spasm_malloc(n * sizeof(int));
  xi = spasm_malloc(n * sizeof(int));
  /* possible optimization : recycle pstack and xi/mark for this */
  rcopy = spasm_malloc((n + 1) * sizeof(int));

  /* first pass */
  for (i = 0; i < n; i++) {
    mark[i] = 0;
  }

  top = n;
  for (i = 0; i < n; i++) {
    if (mark[i]) {
      continue;
    }
    top = spasm_dfs(i, A, top, xi, pstack, mark, SPASM_IDENTITY_PERMUTATION);
  }

  /* second pass : in the order given by xi, in A_t */
  for (i = 0; i < n; i++) {
    mark[i] = 0;
  }

  A_t = spasm_transpose(A, SPASM_IGNORE_VALUES);
  top = n;
  n_scc = n;

  for (k = 0; k < n; k++) {
    i = xi[k];

    if (mark[i]) {
      continue;
    }
    rr[n_scc] = top;
    n_scc--;
    top = spasm_dfs(i, A_t, top, p, pstack, mark, SPASM_IDENTITY_PERMUTATION);
  }

  /* the blocks are now in rr[n_scc:n]. Move them to rr[0:n_scc] */
  rr[n_scc] = 0;
  for (k = n_scc; k <= n; k++) {
    rr[k - n_scc] = rr[k];
    cc[k - n_scc] = rr[k];
  }
  n_scc = n - n_scc;
  P->nr = n_scc;
  P->nc = n_scc;

  /* at this point, blocks are in reverse order in both r and p */

  /* mark rows with (reverse) block numbers */
  for (k = 0; k < n_scc; k++) {
    for (i = rr[k]; i < rr[k + 1]; i++) {
      mark[p[i]] = n_scc - k - 1;
    }
  }

  /* rcopy[k] indicates the next entry in block k in p */
  for (k = 0; k <= n_scc; k++) {
    rcopy[k] = n - rr[n_scc - k];
  }

  /* rcopy[k] indicates the next entry in block k in p */
  for (k = 0; k <= n_scc; k++) {
    rr[k] = rcopy[k];
    cc[k] = rcopy[k];
  }

  for (i = 0; i < n; i++) {
    k = mark[i];                /* the block id */
    p[rcopy[k]] = i;
    q[rcopy[k]] = i;
    rcopy[k]++;
  }

  free(xi);
  free(pstack);
  free(mark);
  free(rcopy);
  spasm_csr_free(A_t);
  spasm_partition_tighten(P);
  return P;
}
