#include <assert.h>
#include "spasm.h"


/* k indicates the starting row.
   imatch[j] indicates which row matches column j (or -1 if column j is not matched).
 */
static int spasm_augmenting_path(const spasm *A, int k, int *row_stack, int *col_stack, int *pointer_stack, int *marks, int *cheap, int *imatch, int *jmatch) {
  int i, j, p, head, found;
  int *Ap, *Aj;

    /* check inputs */
    Ap = A->p;
    Aj = A->j;

    /* initialize the recursion stack (row nodes waiting to be traversed). The
     * stack is held at the begining of xi, and has head elements. */
    head = 0;
    row_stack[head] = k;
    found = 0;

    /* stack empty ? */
    while (head >= 0) {
        /* get i from the top of the recursion stack */
        i = row_stack[head];

	if (marks[i] != k) {
	  marks[i] = k;

	  /* now try to find an unmatched column adjacent to row i */
	  for (p = cheap[i]; p < Ap[i + 1]; p++) {
	    j = Aj[p];
	    if (imatch[j] == -1) {
	      found = 1;
	      break;
	    }
	  }
	  cheap[i] = p;

	  if (found) {
	    col_stack[head] = j;
	    break;
	  }

	  /* all adjacent columns are matched: we need to push the DFS search deeper */
	  pointer_stack[head] = Ap[i];
	}

	/* --- Depth-first-search of column adjacent to row i -------------------- */
	for (p = pointer_stack[head]; p < Ap[i + 1]; p++) {
	  j = Aj[p];
	  if (marks[ imatch[j] ] == k) {
	    continue;
	  }

	  /* found an explored neighbor imatch[j]: pause dfs of node i */
	  pointer_stack[head] = p + 1;
	  col_stack[head] = j;

	  /* start dfs at row imatch[j] */
	  head++;
	  row_stack[head] = imatch[j];
	  break;
        }
	/* row node i is done: pop it from stack */
        if (p == Ap[i + 1]) {
	  head--;
	}
    }

    if (found) {
      //      printf("path of length %d\n", head + 1);
      for (p = head; p >= 0; p--) {
	imatch[ col_stack[p] ] = row_stack[p];
	jmatch[ row_stack[p] ] = col_stack[p];
      }
      return 1;
    }
    return 0;
}


/* returns the size of the matching */
int spasm_maximum_matching(const spasm *A, int *imatch, int *jmatch) {
  int n, m, r, i, j, k, p;
  int *Ap, *Aj, *w, *row_stack, *col_stack, *marks, *pointer_stack, *cheap;

  n = A->n;
  m = A->m;
  r = spasm_min(n, m);
  Ap = A->p;
  Aj = A->j;

  // cette procédure nécessite peut-être que n >= m (?)
  for(j = 0; j < m; j++) {
    imatch[j] = -1;
  }

  /* get workspace */
  w = spasm_malloc(5 * n * sizeof(int));
  row_stack = w;
  col_stack = w + n;
  pointer_stack = w + 2 * n;
  marks = w + 3 * n;
  cheap = w + 4 * n;

  for(i = 0; i < n; i++) {
    marks[i] = -1;
    cheap[i] = Ap[i];
    jmatch[i] = -1;
  }

  k = 0;
  /* --- greedily finds a "cheap" matching -------------------- */
  for(i = 0; i < n; i++) {
    for (p = Ap[i]; p < Ap[i + 1]; p++) {
      j = Aj[p];
      if (imatch[j] == -1) {
	imatch[j] = i;
	jmatch[i] = j;
	k++;
	break;
      }
    }
    cheap[i] = p;
    }

  /* --- finds a maximal matching ----------------------------- */
  for(i = 0; i < n; i++) {
    if (k == r) {
      break; // early abort
    }
    if (jmatch[i] == -1) {
      k += spasm_augmenting_path(A, i, row_stack, col_stack, pointer_stack, marks, cheap, imatch, jmatch);
    }
  }

  free(w);
  return k;
}


int spasm_structural_rank(const spasm *A) {
  int n, m;
  int *imatch, *jmatch;

  n = A->n;
  m = A->m;

  imatch = spasm_malloc(m * sizeof(int));
  jmatch = spasm_malloc(n * sizeof(int));

  return spasm_maximum_matching(A, imatch, jmatch);
}
