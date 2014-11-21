#include <assert.h>
#include "spasm.h"


static inline int spasm_flip(int i) {
    return -i - 2;
}

static inline int spasm_unflip(int i) {
    return (i < 0) ? spasm_flip(i) : i;
}

static inline int spasm_is_marked(const int *w, int j) {
    return (w[j] < 0);
}

static inline void spasm_mark(int *w, int j) {
    w[j] = spasm_flip(w[j]);
}


/*
 * depth-first-search of the graph of a matrix, starting at node j. All nodes
 * encountered during the graph traversal are marked, and added to xi.
 * 
 * j : root of the search G : the graph to search xi : size n. Used both as
 * workspace and to return the result. At the end, the list of traversed
 * nodes is in xi[top:n] pstack : size-n workspace. Used to count the
 * neighbors already traversed.
 * 
 * return value : top
 */
int spasm_dfs(int j, spasm * G, int top, int *xi, int *pstack, const int *pinv) {
    int i, p, p2, done, jnew, head, *Gp, *Gi;

    /* check inputs */
    assert(spasm_is_csc(G));
    assert(xi != NULL);
    assert(pstack != NULL);

    Gp = G->p;
    Gi = G->i;
    /*
     * initialize the recursion stack (nodes waiting to be traversed). The
     * stack is held at the begining of xi, and has head elements.
     */
    head = 0;
    xi[head] = j;

    /* stack empty ? */
    while (head >= 0) {
        /* get j from the top of the recursion stack */
        j = xi[head];
        jnew = (pinv != NULL) ? pinv[j] : j;

        /* has the new node been seen before ? */
        /* neighbors are Gi[Gp[jnew]     : Gp[jnew + 1]] */
        /* UNSEEN neighbors are Gi[pstack[head] : Gp[jnew + 1]] */

        if (!spasm_is_marked(Gp, j)) {
            /* mark node j as seen. This is done only once. */
            spasm_mark(Gp, j);
            /*
             * Initialize pstack for this node: first unseen neighbor is...
             * the first neighbor
             */
            pstack[head] = (jnew < 0) ? 0 : spasm_unflip(Gp[jnew]);
	    /* Note: unflip does nothing if not flipped.
	       So it can be used to simply get the "actual", unaltered value.*/
                }

                /* node j done if no unvisited neighbors */
                done = 1;

            /* index of last neighbor */
            p2 = (jnew < 0) ? 0 : spasm_unflip(Gp[jnew + 1]);
    //Note:same unflip trick

            /* examine all yet-unseen neighbors of j */
                for (p = pstack[head]; p < p2; p++) {

                /* consider next neighbor node, namely i */
                i = Gi[p];
                /* if already visisted, skip */
                if (!spasm_is_marked(Gp, i)) {
                    /*
                     * interrupt the enumeration of neighbors of node jnew,
                     * and deal with i instead.
                     */

                    /* Save number of examined neighbors of jnew */
                    pstack[head] = p;

                    /*
                     * push node i onto the recursion stack. This will start
                     * a DFS from i
                     */
                    ++head;
                    xi[head] = i;
                    /* node j is not done, and exit the loop */
                    done = 0;
                    break;
                }
            }

            /* depth-first search at node j done ? */
            if (done) {
                /* pop j from the recursion stack */
                head--;
                /* and push j in the output stack */
                --top;
                xi[top] = j;
            }
        }
        return top;
    }


    /*
     * Compute the set of nodes from G reachable from a node of B[:,k] (used
     * to determine the pattern of a sparse triangular solve)
     *
     * G : graph to search (modified, then restored)
     * B : RHS (starting point of the search) 
     *k : k-th column of B is used. 
     * xi: size 2n. Used as workspace. Output in xi[top:n] 
     * pinv: mapping of rows to columns of G.
     *
     * return value : top
     *
     * xi [top...n-1] = nodes reachable from graph of G*P'
     *  via nodes in B(:,k).
     *
     * xi [n...2n-1] used as workspace
     */
    int spasm_reach(spasm * G, const spasm * B, int k, int *xi, const int *pinv){
        int p, n, top, *Bp, *Bi, *Gp;

        /* check inputs */
            assert(spasm_is_csc(G));
            assert(spasm_is_csc(B));
            assert(xi != NULL);

            n = G->n;
            Bp = B->p;
            Bi = B->i;
            Gp = G->p;
            top = n;
        /*
         * iterates over the k-th column of B.  For each row index i present
         * in B[:,k], check if i is in the pattern (i.e. if it is marked). If
         * not, start a DFS from i to add all nodes reachable therefrom to
         * the pattern.
         */
        for (p = Bp[k]; p < Bp[k + 1]; p++) {
            if (!spasm_is_marked(Gp, Bi[p])) {
                top = spasm_dfs(Bi[p], G, top, xi, xi + n, pinv);
            }
        }
        /* restore G : unmark all marked nodes. */
        for (p = top; p < n; p++) {
            spasm_mark(Gp, xi[p]);
        }
        return top;
    }
