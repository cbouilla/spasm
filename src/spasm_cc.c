#include <assert.h>
#include "spasm.h"

/* returns the number of non-trivial (size > 0) connected component of A, seen as an
   undirected bipartite graph.

   if jmatch == NULL, columns are sorted in natural order.

   if jmatch != NULL, all rows must be matched. matched columns come
   first (in the order given by the matching), unmatched columns come
   next (in natural order).

   if imatch == NULL, rows are sorted in natural order.

   if imatch != NULL, all columns must be matched. matched rows come
   first (in the order given by the matching), unmatched rows come
   next (in natural order).

   both jmatch != NULL and imatch != NULL makes no sense.

   If the transpose of A is not given, it will be computed.
*/
spasm_partition * spasm_connected_components(const spasm *A, const spasm *givenA_t, const int *jmatch, const int *imatch) {
  int n, m, i, j, k, root, n_cc, rhead, rtail, chead, ctail, px;
  int *Ap, *Aj, *A_tp, *A_tj, *rmark, *cmark, *p, *q, *rr, *cc, *rcopy, *ccopy;
  spasm *A_t;
  spasm_partition *P;

  assert(A != NULL);
  assert(jmatch == NULL || imatch == NULL);

  n = A->n;
  m = A->m;
  A_t = (givenA_t != NULL) ? (spasm *) givenA_t : spasm_transpose(A, SPASM_IGNORE_VALUES);
  Ap = A->p;
  Aj = A->j;
  A_tp = A_t->p;
  A_tj = A_t->j;

  rmark = spasm_malloc(n * sizeof(int));
  cmark = spasm_malloc(m * sizeof(int));
  spasm_init_vector(rmark, n, -1);
  spasm_init_vector(cmark, m, -1);

  P = spasm_partition_alloc(n, m, n + 1, n + 1);
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

    // skip row if it is empty
    if (Ap[root] == Ap[root + 1]) {
      continue;
    }

    /* previous block stops here */
    rr[n_cc] = rhead;
    cc[n_cc] = chead;

    /* start BFS from row root */
    p[rtail] = root;
    rmark[root] = n_cc;
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

	cmark[j] = n_cc;
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
	    continue;
	  }

	  rmark[i] = n_cc;
	  p[rtail] = i;
	  rtail++;
	}
      }
    }

    n_cc++;
  }

  rr[n_cc] = rhead; // not necessarily n
  cc[n_cc] = chead; // not necessarily m (cf. below)

  /* WARNING : DO **NOT** CREATE EMPTY CC */

  P->nr = n_cc;
  P->nc = n_cc;

  /* rows and column are in a somewhat random order in p and q.  put
     them in the natural order inside each block. */

  rcopy = spasm_malloc((n_cc + 1) * sizeof(int));
  ccopy = spasm_malloc((n_cc + 1) * sizeof(int));

  /* rcopy[k] (resp. ccopy[k]) indicates the next row (resp. column) of block k should land in p (resp. q) */
  for(i = 0; i <= n_cc; i++) {
    rcopy[i] = rr[i];
    ccopy[i] = cc[i];
  }

  /* first pass, only unmatched (non-empty) rows (if jmatch != NULL, of course) */
  for(i = 0; i < n; i++) {
    k = rmark[i];
    if (k == -1 || (jmatch != NULL && jmatch[i] != -1)) {
      continue;
    }
    p[ rcopy[k] ] = i;
    rcopy[k]++;
  }

  /* second pass, only matched rows */
  if (jmatch != NULL) {
    for(i = 0; i < n; i++) {
      if (jmatch[i] == -1) {
	continue;
      }
      k = rmark[i];
      p[ rcopy[k] ] = i;
      rcopy[k]++;

      j = jmatch[i];
      cmark[j] = -2; // do not add this column a second time later on
      q[ ccopy[k] ] = j;
      ccopy[k]++;
    }
  }

  /* add empty rows in natural order */
  for(i = 0; i < n; i++) {
    k = rmark[i];
    if (k != -1) {
      continue;
    }
    p[ rcopy[n_cc] ] = i;
    rcopy[n_cc]++;
  }

  /* add remaining non-empty column in natural order */
  for(j = 0; j < m; j++) {
    k = cmark[j];
    if (k < 0) {
      continue;
    }
    q[ ccopy[k] ] = j;
    ccopy[k]++;
  }

  /* add remaining empty column in natural order */
  for(j = 0; j < m; j++) {
    k = cmark[j];
    if (k != -1) {
      continue;
    }
    q[ ccopy[n_cc] ] = j;
    ccopy[n_cc]++;
  }


  /* cleanup */
  free(rcopy);
  free(ccopy);
  free(rmark);
  free(cmark);
  if (givenA_t == NULL) {
    spasm_csr_free(A_t);
  }
  spasm_partition_tighten(P);
  return P;
}
