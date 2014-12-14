#include <assert.h>
#include "spasm.h"

/* returns the number of non-trivial connected component of A, seen as an
   undirected bipartite graph. */
int spasm_connected_components(const spasm *A, int *rmark, int *cmark) {
  int n, m, i, j, p, root, size, n_cc, rhead, rtail, chead, ctail;
  int *Ap, *Aj, *A_tp, *A_tj, *rqueue, *cqueue;
  spasm *A_t;

  n = A->n;
  m = A->m;
  A_t = spasm_transpose(A, 0);
  Ap = A->p;
  Aj = A->j;
  A_tp = A_t->p;
  A_tj = A_t->j;

  rqueue = spasm_malloc(n * sizeof(int));
  cqueue = spasm_malloc(m * sizeof(int));

  for (i = 0 ; i < n ; i++) {
    rmark[i] = -1;
  }
  for (j = 0 ; j < m ; j++) {
    cmark[j] = -1;
  }

  n_cc = 0;
  for(root = 0; root < n; root++) {
    /* identify connected component containing row i, if not already done */

    if (rmark[root] != -1) {
      assert(rmark[root] < root);
      continue;
    }

    /* start BFS */
    rhead = 0;
    rqueue[0] = root;
    rmark[root] = root;
    size = 1;
    rtail = 1;
    chead = 0;
    ctail = 0;

    /* while row queue is not empty */
    while (rhead < rtail) {
      i = rqueue[rhead];
      rhead++;

      assert(rmark[i] == root);

      for (p = Ap[i]; p < Ap[i + 1]; p++) {
	j = Aj[p];
	if (cmark[j] != -1) {
	  assert(cmark[j] == root);
	  continue;
	}

	cmark[j] = root;
	cqueue[ctail] = j;
	ctail++;
      }

      /* while col queue is not empty */
      while (chead < ctail) {
	j = cqueue[chead];
	chead++;

	for (p = A_tp[j]; p < A_tp[j + 1]; p++) {
	  i = A_tj[p];
	  if (rmark[i] != -1) {
	    assert(rmark[i] == root);
	    continue;
	  }

	  rmark[i] = root;
	  rqueue[rtail] = i;
	  rtail++;
	  size++;
	}
      }
    }

    n_cc += (size > 1);
  }

  free(rqueue);
  free(cqueue);
  spasm_csr_free(A_t);
  return n_cc;
}
