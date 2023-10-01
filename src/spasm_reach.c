#include <assert.h>

#include "spasm.h"

/**
 * Depth-first-search along alternating paths of a bipartite graph.
 *
 * If a column j is pivotal (qinv[j] != -1), then move to the row (call it i)
 * containing the pivot; explore columns adjacent to row i, depth-first. 
 * The traversal starts at column jstart.
 *
 * qinv[j] indicates the row on which the j-th column pivot can be found.
 * qinv[j] == -1 means that there is no pivot on column j.
 *
 * xj is of size m (#columns). Used both as workspace and to return the result.
 * At the end, the list of traversed nodes is in xj[top:m].  This returns top.
 *
 * pstack : size-m workspace (used to count the neighbors already traversed)
 * marks  : size-m workspace (indicates which columns have been seen already)
 */
int spasm_dfs(int jstart, const struct spasm_csr *A, int top, int *xj, int *pstack, int *marks, const int *qinv)
{
	/* check inputs */
	assert(A != NULL);
	assert(xj != NULL);
	assert(pstack != NULL);
	assert(marks != NULL);
	assert(qinv != NULL);

	const i64 *Ap = A->p;
	const int *Aj = A->j;
	/*
	 * initialize the recursion stack (columns waiting to be traversed).
	 * The stack is held at the begining of xj, and has head elements.
	 */
	int head = 0;
	xj[head] = jstart;

	/* stack empty ? */
	while (head >= 0) {
		/* get j from the top of the recursion stack */
		int j = xj[head];
		int i = qinv[j];       /* row with the pivot on column j, or -1 if none */

		if (!marks[j]) {
			/* mark column j as seen and initialize pstack. This is done only once. */
			marks[j] = 1;
			pstack[head] = 0;
		}

		if (i < 0) {
			/* push initial column in the output stack and pop it from the recursion stack*/
			top -= 1;
			xj[top] = xj[head];
			head -= 1;
			continue;
		}

		/* size of row i */
		int p2 = spasm_row_weight(A, i);

		/* examine all yet-unseen entries of row i */
		int k;
		for (k = pstack[head]; k < p2; k++) {
			i64 px = Ap[i] + k;
			int j = Aj[px];
			if (marks[j])
				continue;
			/* interrupt the enumeration of entries of row i, and start DFS from column j */
			pstack[head] = k + 1;   /* Save status of row i in the stack. */
			xj[++head] = j;         /* push column j onto the recursion stack */
			break;
		}
		if (k == p2) {
			/* row i fully examined; push initial column in the output stack and pop it from the recursion stack */
			top -= 1;
			xj[top] = xj[head];
			head -= 1;
		}
	}
	return top;
}


/*
 * Reachability along alternating paths of a bipartite graph.
 * Compute the set of columns of A reachable from all columns indices in B[k]
 * (this is used to determine the pattern of a sparse triangular solve)
 * 
 * xj must be preallocated of size 3*m and zeroed out on entry.
 * On output, the set of reachable columns is in xj[top:m].
 * This returns top.  xj remains in a usable state (no need to zero it out again)
 *
 * qinv locates the pivots in A.
 *
 * This function does not require the pivots to be the first entries of the rows.
 */
int spasm_reach(const struct spasm_csr *A, const struct spasm_csr *B, int k, int l, int *xj, const int *qinv)
{
	/* check inputs */
	assert(A != NULL);
	assert(B != NULL);
	assert(xj != NULL);
	assert(qinv != NULL);

	const i64 *Bp = B->p;
	const int *Bj = B->j;
	int m = A->m;
	int top = m;
	int *pstack = xj + m;
	int *marks = pstack + m;

	/*
	 * iterates over the k-th row of B.  For each column index j present
	 * in B[k], check if j is in the pattern (i.e. if it is marked). If
	 * not, start a DFS from j and add to the pattern all columns
	 * reachable from j.
	 */
	for (i64 px = Bp[k]; px < Bp[k + 1]; px++) {
		int j = Bj[px];
		if (!marks[j])
			top = spasm_dfs(j, A, top, xj, pstack, marks, qinv);
	}

	/* unmark all marked nodes. */
	/*
	 * TODO : possible optimization : if stuff is marked "k", and
	 * initialized with -1, then this is not necessary
	 */
	for (int px = top; px < l; px++) {
		int j = xj[px];
		marks[j] = 0;
	}
	return top;
}
