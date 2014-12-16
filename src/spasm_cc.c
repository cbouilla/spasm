#include <assert.h>
#include "spasm.h"

/* returns the number of non-trivial (size > 0) connected component of A, seen as an
   undirected bipartite graph. */
spasm_partition * spasm_connected_components(const spasm *A) {
  int n, m, i, j, root, n_cc, rhead, rtail, chead, ctail, px;
  int *Ap, *Aj, *A_tp, *A_tj, *rmark, *cmark, *p, *q, *rr, *cc;
  spasm *A_t;
  spasm_partition *P;

  n = A->n;
  m = A->m;
  A_t = spasm_transpose(A, 0);
  Ap = A->p;
  Aj = A->j;
  A_tp = A_t->p;
  A_tj = A_t->j;

  rmark = spasm_malloc(n * sizeof(int));
  cmark = spasm_malloc(m * sizeof(int));
  for (i = 0 ; i < n ; i++) {
    rmark[i] = -1;
  }
  for (j = 0 ; j < m ; j++) {
    cmark[j] = -1;
  }

  P = spasm_partition_alloc(n, m, n, m);
  p = P->p;
  q = P->q;
  rr = P->rr;
  cc = P->cc;
  rhead = 0;
  rtail = 0;
  chead = 0;
  ctail = 0;
  n_cc = 0;

  for(root = 0; root < n; root++) {

    if (rmark[root] != -1) {
      continue;
    }

    /* previous block stops here */
    rr[n_cc] = chead;
    cc[n_cc] = chead;

    /* start BFS from row root */
    p[rtail] = root;
    rmark[root] = root;
    rtail++;

    /* while row queue is not empty */
    while (rhead < rtail) {
      i = p[rhead];
      rhead++;

      for (px = Ap[i]; px < Ap[i + 1]; px++) {
	j = Aj[px];
	if (cmark[j] != -1) {
	  continue;
	}

	cmark[j] = root;
	q[ctail] = j;
	ctail++;
      }

      /* while col queue is not empty */
      while (chead < ctail) {
	j = q[chead];
	chead++;

	for (px = A_tp[j]; px < A_tp[j + 1]; px++) {
	  i = A_tj[px];
	  if (rmark[i] != -1) {
	    assert(rmark[i] == root);
	    continue;
	  }

	  rmark[i] = root;
	  p[rtail] = i;
	  rtail++;
	}
      }
    }

    n_cc++;
  }

  rr[n_cc] = n;
  cc[n_cc] = m;
  P->nr = n_cc;
  P->nc = n_cc;

  free(rmark);
  free(cmark);
  spasm_csr_free(A_t);

  spasm_partition_tighten(P);
  return P;
}
