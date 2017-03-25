#include <assert.h>
#include "spasm.h"


/*given a matrix A, a row i of A, a positive integer j (the father of i) such that Aij != 0, 
 *search j in row i. When j is found swap j and Aj[At[i]], At[i] <- At[i]+1.
 *.
 */
void spasm_search_father(spasm *A, int i, int j, int *At){
  int *Ap, *Aj, *Ax, px;

  assert(A != NULL);
  assert(At != NULL);
  Ap = A->p;
  Aj = A->j;
  Ax = (A->x != NULL) ? A->x: NULL;

  for(px = At[i]; px < Ap[i+1]; px++){ // j isn't in the tree yet, so it is after At[i].
    if(Aj[px] == j){
      break; // once j is found, quit the loop.
    }
  }

  assert(Aj[px] == j); // Aij != 0.
  spasm_swap(Aj, px, At[i]);
  if(A->x != NULL){
    spasm_swap(Ax, px, At[i]);
  }

  //update At.
  At[i]++;

  return;
}

/*given a table match over rows and the cols, a row i and a col j, check wether i and j are in the
 *matching. if not mark (i,j) as an edge of the matching.
 *
 *i: index of a row
 *j: index of a column
 *
 *match: table of the matching (match.r[i] = j; and match.c[j] = i if (ij) is in the matching).
 *
 *return 1 if (i,j) is a valid vertice of the matching and 0 if not.
 */
int spasm_tree_update_matching(spasm_rc *match, int i, int j, int b){

  if((match->r[i] == -1) && (match->c[j] == -1)){
    match->r[i] = j;
    match->c[j] = i;
    return 1;
  }

  return 0;
}

/*
 *bipartite dfs:
 *A: input matrix.
 *tA: transpose of input matrix (not as arg)
 *root: root of the tree (a row index)
 *first_passage (fpr and fpc): init to -1, 0, when first passage, mark the number of backward edge. 
 *match: store the matching
 *At: mark the end of the index that are in the tree in each row/col.
 *
 *return 1 if all nodes have been watch at least once, 0 otherwise.
 */
