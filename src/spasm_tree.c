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

  assert(At[i] == Ap[i]); // first passage on this node.

  for(px = At[i]; px < Ap[i+1]; px++){ 
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



/*
 *bipartite dfs:
 *A: input matrix.
 *tA: transpose of input matrix
 *root: root of the tree (a row index)
 *first_passage (fpr and fpc): init to -1, 0, when first passage, mark the number of backward edge. 
 *match: store the matching
 *At: mark the end of the index that are in the tree in each row/col.
 *
 *return 1 if all nodes have been watch at least once, 0 otherwise.
 */
int spasm_dfs_bipartite(spasm *A, spasm *TA, int root, spasm_rc *first_passage, spasm_rc *match, spasm_rc *At){
  int *Ap, *tAp, *Aj, *tAj, *pstack, *Atr, *Atc, *matchr, *matchc, *fpr, *fpc;
  int j, head, father, px, n, m, count;

  assert(A != NULL);
  assert(TA != NULL);
  assert(first_passage != NULL);
  assert(match != NULL);
  assert(At != NULL);
  Ap = A->p;
  Aj = A->j;
  tAp = TA->p;
  tAj = TA->j;
  Atr = At->r;
  Atc = At->c;
  matchr = match->r;
  matchc = match->c;
  fpr = first_passage->r;
  fpc = first_passage->c;
  n = A->n;
  m = A->m;
  assert(n == TA->m);
  assert(m == TA->n);

  //init variables.
  head = 0;
  father = -1; // root has no father.
  count = 0;

  //get workspace
  pstack = spasm_malloc((n+m)*sizeof(int));

  //init workspace
  pstack[head] = root; // stack root. root is a row, so it is >0.

  //main loop: dfs
  while(head >= 0){

    if(pstack[head] <0){ // pstack[head] is a column

      assert(father >= 0); // root isn't a column, so there is a father.

      j = - pstack[head] - 1; // get the id of the column.
      if(fpc[j] == -1){ // j has not been seen yet.

	count++; // a new node is seen.

	//is (father, j) an edge of the matching?
	assert(matchc[j] == -1); // j has not been seen yet.
	if(matchr[father] == -1){ // if father isn't a node of the matching.
	  matchr[father] = j;
	  matchc[j] = father; // update matching.
	}

	spasm_search_father(TA, j, father, Atc); //put father at the end of the tree;
	fpc[j]++; //mark j.

      }else{ // j is marked and no backward edge possible => change of branch.

	father = tAj[tAp[j]]; // j is a column so the first entry of tAp is the father.
      }

      // search for descendance.
      for(px = Atc[j]; px < tAp[j+1]; px++){
	
	if(fpr[tAj[px]] != -1){ // row i = tAj[px] is already in the tree.
	  fpr[tAj[px]] ++; // new backward edge.
	  continue;
	}
	
	father = j; // j is the father of tAj[px].
	j = tAj[px];
	
	//put the new node at the end of the tree.
	spasm_swap(tAj, px, Atc[j]); // no numerical value needed for transpose.
	Atc[j]++; // new node in the tree.
	
	//update stack.
	head++;
	pstack[head] = j; // j is a row.
	break;
	
      }

	//stop test: end of the branch
	if(px == tAp[j+1]){
	  //j has no son.
	  j =  father; // go back in the tree.
	  head--;
	}
	
    } else{ // pstack[head] is a row.

      j = pstack[head];

      if((fpr[j] == -1) && (j==root)){ // first passage on root.

	count++; // new node.
	
	assert(father == -1);
	fpr[j]++; // mark root.
	
      } else if((fpr[j] == -1) && (j != root)){ // j has not been seen yet.
	//is (father, j) an edge of the matching?
	assert(father >= 0); // j is not root.

	count++;
	
	assert(matchr[j] == -1); // j has not been seen yet.
	if(matchc[father] == -1){ // if father isn't a node of the matching.
	  matchc[father] = j;
	  matchr[j] = father; // update matching.
	}

	spasm_search_father(A, j, father, Atr); //put father at the begining of the row;
	fpr[j]++; //mark j.

      }else if(j != root){ // j is marked and no backward edge possible => change of branch.

	father = Aj[Ap[j]]; // j is not root so the first entry of Ap is the father.
	
      }else{ // j is root and is marked.

	father = -1; // no father.
      }

      // search for descendance.
      for(px = Atr[j]; px < Ap[j+1]; px++){
	
	if(fpc[Aj[px]] != -1){ // col i = Aj[px] is already in the tree.
	  fpr[Aj[px]] ++; // new backward edge.
	  continue;
	}

	//One son found:
	
	father = j; // j is the father of Aj[px].
	j = Aj[px];
	
	//put the new node at the end of the tree.
	spasm_swap(Aj, px, Atr[j]);
	if(A->x != NULL){
	  spasm_swap(A->x, px, Atr[j]); // don't forget numerical values.
	}
	Atr[j]++;
	
	//update stack.
	head++;
	pstack[head] = j; // j is a row.
	break;
	
      }

	//stop test: end of the branch
	if(px == tAp[j+1]){
	  //j has no son.
	  j =  father; // go back in the tree.
	  head--;
	}
      
    
    }
  }

  //free workspace
  free(pstack);

  return ((count) == (n+m));
}
