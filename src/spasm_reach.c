#include <assert.h>
#include "spasm.h"

/**
 * depth-first-search of the graph of a matrix, starting at column j. 
 *
 * The traversal works as follow : starting from a column j, find the
 * row containing the pivot on this column (row "inew"). If there is no
 * pivot on column j, then assume it is on row j (this allows to work
 * with matrices like U for which n < m --- an implicit identity is padded below).
 *
 * When pivot row has been found, append its next unseen entry to the stack, and restart.
 *
 * @param pinv 
 *           pinv[j] indicates the row on which the j-th column pivot can be found.
 *           pinv[j] == -1 means that there is no pivot on column j.
 *      
 * The matrix need not be square. If n < m, then pinv must be non-NULL, and
 * pinv[i] = -1 means that row[i] implicitly is an identity row.
 *
 * @param j column at which the traversal starts.
 *
 * @param G the graph to search
 *
 * @param xj : size m. Used both as workspace and to return the result. At the end, the
 * list of traversed nodes is in xj[top:m]
 *
 * @param pstack 	size-n workspace. Used to count the neighbors already traversed.
 *
 * @marks : size-m workspace. Indicates which column nodes have been dealt with.
 *
 * @return top
 */
int spasm_dfs(int j, const spasm * G, int top, int *xj, int *pstack, int *marks, const int *pinv) {
    int p, p2, inew, head, *Gp, *Gj;

    /* check inputs */
    assert(G != NULL);
    assert(xj != NULL);
    assert(pstack != NULL);
    assert(marks != NULL);

    Gp = G->p;
    Gj = G->j;
    /*
     * initialize the recursion stack (rows waiting to be traversed). The
     * stack is held at the begining of xj, and has head elements.
     */
    head = 0;
    xj[head] = j;

    /* stack empty ? */
    while (head >= 0) {
	/* get i from the top of the recursion stack */
	j = xj[head];
	inew = (pinv != NULL) ? pinv[j] : j;

	/*
	 * has row i been seen before ?
	 * adjacent columns are Gj[ Gp[jnew]     : Gp[jnew + 1] ]
	 * UNSEEN columns are   Gj[ pstack[head] : Gp[jnew + 1] ]
	 */

	if (!marks[j]) {
	    /* mark node i as seen. This is done only once. */
	    marks[j] = 1;
	    /*
	     * Initialize pstack for this node: first unseen column is... the
	     * first entry on the row
	     */
	    pstack[head] = (inew < 0) ? 0 : Gp[inew];
	}
	/* index of last entry */
	p2 = (inew < 0) ? 0 : Gp[inew + 1];

	/* examine all yet-unseen entries of row i */
	for (p = pstack[head]; p < p2; p++) {

	    /* consider next adjacent column */
	    j = Gj[p];

	    /* if already visisted, skip */
	    if (marks[j]) {
		continue;
	    }
	    /*
	     * interrupt the enumeration of entries of row inew, and deal
	     * with column j instead. Save status of row inew in the stack.
	     */
	    pstack[head] = p + 1;

	    /*
	     * push column j onto the recursion stack. This will start a DFS
	     * from j
	     */

	    head++;
	    xj[head] = j;

	    /* node i is not done, exit the loop */
	    break;
	}

	/* depth-first search at node i done ? */
	if (p == p2) {
	    /* push initial column in the output stack */
	    top--;
	    xj[top] = xj[head]; 
	    
	    /* pop it from the recursion stack */
	    head--;
	}
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
int spasm_reach(const spasm * G, const spasm * B, int k, int l, int *xj, const int *pinv) {
    int p, m, top, *Bp, *Bj, *pstack, *marks;

    /* check inputs */
    assert(G != NULL);
    assert(B != NULL);
    assert(xj != NULL);

    m = G->m;
    Bp = B->p;
    Bj = B->j;
    top = l;

    pstack = xj + l;
    marks = pstack + l;

    /*
     * iterates over the k-th row of B.  For each column index j
     * present in B[k], check if j is in the pattern (i.e. if it is
     * marked). If not, start a DFS from j and add to the pattern all
     * columns reachable from j.
     */
    for (p = Bp[k]; p < Bp[k + 1]; p++) {
	if (!marks[Bj[p]]) {
	    top = spasm_dfs(Bj[p], G, top, xj, pstack, marks, pinv);
	}
    }

    /* unmark all marked nodes. */
    for (p = top; p < l; p++) {
      marks[ xj[p] ] = 0;
    }

    return top;
}


/*
 * Compute the set of nodes from G reachable from any node in B[k] (used to
 * determine the pattern of a sparse triangular solve)
 *
 * G : graph to search
 *
 * yi : rhs pattern.
 *
 * start : begining of yi. end : end of yi.
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
int spasm_scat_reach(const spasm *G, int *yj, int start, int end, int l, int *xj, const int *pinv){
int p, top, *pstack, *marks;

    /* check inputs */
    assert(G != NULL);
    assert(yj != NULL);
    assert(xj != NULL);

    top = l;

    pstack = xj + l;
    marks = pstack + l;

    /*
     * iterates over y.  For each column index j
     * present in yi, check if j is in the pattern (i.e. if it is
     * marked). If not, start a DFS from j and add to the pattern all
     * columns reachable from j.
     */
    
    for (p = start; p < end; p++) {
	if (!marks[yj[p]]) {
	    top = spasm_dfs(yj[p], G, top, xj, pstack, marks, pinv);
	}
    }

    /* unmark all marked nodes. */
    for (p = top; p < l; p++) {
      marks[ xj[p] ] = 0;
    }

    return top;
}
