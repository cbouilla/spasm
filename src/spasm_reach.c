#include <assert.h>
#include "spasm.h"

/*
 * depth-first-search of the graph of a matrix, starting at node j. All nodes
 * encountered during the graph traversal are marked, and added to xi.
 *
 * The matrix need not be square. If n < m, then pinv must be non-NULL,
 *   and pinv[i] = -1 means that row[i] implicitly is an identity row.
 *
 * i : root of the search
 *
 * G : the graph to search
 *
 * xi : size m. Used both as workspace and to return the result.
 * At the end, the list of traversed nodes is in xi[top:m]
 *
 * pstack : size-n workspace. Used to count the neighbors already traversed.
 *
 * marks : size-m
 *
 * return value : top
 */
int spasm_dfs(int i, const spasm * G, int top, int *xi, int *pstack, int *marks, const int *pinv) {
    int j, p, p2, inew, head, *Gp, *Gj;

    /* check inputs */
    assert(G != NULL);
    assert(xi != NULL);
    assert(pstack != NULL);
    assert(marks != NULL);

    Gp = G->p;
    Gj = G->j;
    /*
     * initialize the recursion stack (rows waiting to be traversed). The
     * stack is held at the begining of xi, and has head elements.
     */
    head = 0;
    xi[head] = i;

    /* stack empty ? */
    while (head >= 0) {
        /* get i from the top of the recursion stack */
        i = xi[head];
        inew = (pinv != NULL) ? pinv[i] : i;

        /* has row i been seen before ?
         * adjacent columns are Gj[     Gp[jnew] : Gp[jnew + 1] ]
         * UNSEEN columns are   Gj[ pstack[head] : Gp[jnew + 1] ] */

        if (!marks[i]) {
            /* mark node i as seen. This is done only once. */
  	    marks[i] = 1;
            /*
             * Initialize pstack for this node: first unseen column is...
             * the first entry on the row
             */
            pstack[head] = (inew < 0) ? 0 : Gp[inew];
        }

	/* index of last entry */
	p2 = (inew < 0) ? 0 : Gp[inew + 1];

	/* examine all yet-unseen entries of row i */
	for (p = pstack[head]; p < p2; p++) {

	  /* consider next adjacent column, namely j */
	  j = Gj[p];

	  /* if already visisted, skip */
	  if (marks[j]) {
	    continue;
	  }

	  /* interrupt the enumeration of neighbors of node inew,
	     and deal with j instead. Save index of examined neighbors of inew */
	  pstack[head] = p + 1;

	  /* push node j onto the recursion stack. This will start a DFS from j */

	  head++;
	  xi[head] = j;
	  /* node i is not done, and exit the loop */
	  break;
	}

	/* depth-first search at node i done ? */
	if (p == p2) {
	  /* pop i from the recursion stack */
	  head--;
	  /* and push i in the output stack */
	  top--;
	  xi[top] = i;
	}
    }
    return top;
}


    /*
     * Compute the set of nodes from G reachable from any node in B[k] (used
     * to determine the pattern of a sparse triangular solve)
     *
     * G : graph to search (modified, then restored)
     *
     * B : RHS (starting point of the search)
     *
     * k : k-th column of B is used.
     *
     * xi: size 3m. Used as workspace. Output in xi[top:m]
     *
     * pinv: mapping of rows to columns of G.
     *
     * return value : top
     *
     * xi [top...m-1] = nodes reachable from graph of G*P'
     *  via nodes in B(:,k).
     *
     * xi [m...3n-1] used as workspace
     */
    int spasm_reach(const spasm * G, const spasm * B, int k, int *xi, const int *pinv) {
      int p, m, top, *Bp, *Bj, *Gp, *pstack, *marks;

        /* check inputs */
            assert(G != NULL);
            assert(B != NULL);
            assert(xi != NULL);

            m = G->m;
            Bp = B->p;
            Bj = B->j;
            Gp = G->p;
            top = m;

	    pstack = xi + m;
	    marks = pstack + m;
	    for (p = 0; p < m; p++) {
	      marks[p] = 0;
	    }

        /*
         * iterates over the k-th row of B.  For each column index j present
         * in B[k], check if i is in the pattern (i.e. if it is marked). If
         * not, start a DFS from i to add all nodes reachable therefrom to
         * the pattern.
         */
        for (p = Bp[k]; p < Bp[k + 1]; p++) {
            if (!marks[ Bj[p] ]) {
	      top = spasm_dfs(Bj[p], G, top, xi, pstack, marks, pinv);
            }
        }
        /* restore G : unmark all marked nodes. */
        return top;
    }
