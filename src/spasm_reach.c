#include <assert.h>

#include "spasm.h"

/**
 * depth-first-search of the graph of a matrix, starting at column j.
 *
 * The traversal works as follow : starting from a column j, find the
 * row containing the pivot on this column (row "inew"). If there is not pivot on column j, stop.
 *
 * When a pivot row has been found, append its next unseen entry to the stack, and restart.
 *
 * @param qinv
 *           qinv[j] indicates the row on which the j-th column pivot can be found.
 *           qinv[j] == -1 means that there is no pivot on column j.
 *
 * The matrix need not be square. If n < m, then qinv must be non-NULL, and
 * qinv[j] = -1 that there is no pivot on column j
 *
 * @param j column at which the traversal starts.
 *
 * @param G the graph to search
 *
 * @param xj : size m. Used both as workspace and to return the result. At the end, the
 * list of traversed nodes is in xj[top:m]
 *
 * @param pstack size-n workspace. Used to count the neighbors already traversed.
 *
 * @marks : size-m workspace. Indicates which column nodes have been dealt with.
 *
 * @return top
 */
int spasm_dfs(int start, const spasm *A, int top, int *xj, int *pstack, int *marks, const int *qinv)
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
	xj[head] = start;

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
			xj[--top] = xj[head--];
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
		if (k == p2)
			/* row i fully examined; push initial column in the output stack and pop it from the recursion stack */
			xj[--top] = xj[head--];
	}
	return top;
}


/*
 * Compute the set of nodes from G reachable from any node in B[k] (used to
 * determine the pattern of a sparse triangular solve)
 * 
 * G : graph to search
 * 
 * B : RHS (starting point of the search)
 * 
 * k : k-th row of B is used.
 * 
 * l : upper-bound on both dimensions of the matrix
 * 
 * xj: size 3l. Used as workspace. Output in xj[top:l]
 * 
 * pinv: mapping of rows to columns of G.
 * 
 * return value : top
 * 
 * xj [top...l-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * 
 * xj [l...3l-1] used as workspace
 */
int spasm_reach(const spasm *A, const spasm *B, int k, int l, int *xj, const int *qinv)
{
	/* check inputs */
	assert(A != NULL);
	assert(B != NULL);
	assert(xj != NULL);
	assert(qinv != NULL);

	const i64 *Bp = B->p;
	const int *Bj = B->j;
	int top = l;
	int *pstack = xj + l;
	int *marks = pstack + l;

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
