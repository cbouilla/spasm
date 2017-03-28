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
 *col_TO: table that sort the columns in topological order.
 *
 *return 1 if all nodes have been watch at least once, 0 otherwise.
 */
int spasm_dfs_bipartite(spasm *A, spasm *TA, int root, spasm_rc *first_passage, spasm_rc *match, spasm_rc *At, int *col_TO, int count, int *count_col_pt){
  int *Ap, *tAp, *Aj, *tAj, *pstack, *Atr, *Atc, *matchr, *matchc, *fpr, *fpc;
  int j, head, father, px, n, m, count_col;


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
  count_col = *count_col_pt;
  
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
	col_TO[count_col] = j;
	count_col++;

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
	spasm_swap(tAj, px, Atc[father]); // no numerical value needed for transpose.
	Atc[father]++; // new node in the tree.
	
	//update stack.
	head++;
	
	pstack[head] = j; // j is a row.
	fprintf(stderr,"pstack[%d] = %d\n", head, pstack[head]);
	break;
	
      }

	//stop test: end of the branch
      if(px == tAp[j+1] ){
	  //j has no son.
	  j =  father; // go back in the tree.
	  head--;
	  fprintf(stderr,"pstack[%d] = %d\n", head, pstack[head]);
	  
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
	  fpr[Aj[px]]++; // new backward edge.
	  continue;
	}

	//One son found:
	
	father = j; // j is the father of Aj[px].
	j = Aj[px];
	
	//put the new node at the end of the tree.
	spasm_swap(Aj, px, Atr[father]);
	if(A->x != NULL){
	  spasm_swap(A->x, px, Atr[father]); // don't forget numerical values.
	}
	Atr[father]++;
	
	//update stack.
	head++;
	pstack[head] = -j -1; // j is a col.
	fprintf(stderr,"pstack[%d] = %d\n", head, pstack[head]);
	break;
	
      }
      
	//stop test: end of the branch
	if(px == Ap[j+1]){
	  //j has no son.
	  j =  father; // go back in the tree.
	  head--;
	  fprintf(stderr,"pstack[%d] = %d\n", head, pstack[head]);
	 
	}
	 
    }
   
  }

  *count_col_pt = count_col;
  //free workspace
  free(pstack);
  //fprintf(stderr, "count vs tot: %d vs %d\n", count, A->n +A->m);
  return count;
}


/*
 * procedure that test whether there is an alterning cycle in O(1).
 */
/* void spasm_matching_mark(spasm *A, spasm *TA, int *cpt, spasm_rc *At, int i, int k, int start, int in_match, int r){ */
/*   int *Atr, *Atc, *Ap, *tAp, *Aj, *tAj; */
/*   int j, px; */
  
/*   assert(A != NULL); */
/*   assert(TA != NULL); */
/*   Atr = At->r; */
/*   Atc = At->c; */
/*   Ap = A->p; */
/*   Aj = A->j; */
/*   tAp = TA->p; */
/*   tAj = TA->j; */
  
/*   cpt[i] = k; // init cpt[u]. */
  
/* } */


/*
 * given two index son and root, go up in the tree from "son" to "root", and check whether there are two
 * successive non matching edge. if so stop the procedure and return 1. if you reach root, return 0.
 */
int spasm_go_up_in_tree(spasm *A, spasm *TA, int son, int root, spasm_rc *match, int r){
  int *Ap, *tAp, *Aj, *tAj, *matchr, *matchc;
  int father, in_match;

  
  assert(A != NULL);
  assert(TA != NULL);

  Ap = A->p;
  tAp = TA->p;
  Aj = A->j;
  tAj = TA->j;
  matchr = match->r;
  matchc = match->c;

  
  father = (r == 1) ? Aj[Ap[son]] : tAj[tAp[son]]; //father is the first entry on row son.
 

  //The first edge isn't in the matching.
  assert(matchr[son] != father);
  in_match = 0;
  
  while(father != root){
    r ^= 1; // <-- if son was a col, son is now a row and resp.
    son = father;
    father = (r == 1) ? Aj[Ap[son]] : tAj[tAp[son]];

    if(in_match == 0){ // last edge tested isn't in the matching.
      if(((r == 1) && (matchr[son] == father)) || ((r == 0) && matchc[son] == father)){
	in_match = 1;
      }else{ // two successive edge not in the matching.
	return 1; 
      }
    }
    
  }
  
  return 0; // root has been reached. alterning path exist.
}

/*
 * Given u such that (u,v) is in the matching. Test whether (u,v) is valid (no alterning cycle).
 * Note, we only need to test whether the father u has backward edge that goes to the bottom.
 *
 * Test all "critical" nodes (i.e. nodes w, st (u,w) is backward). 
 *
 * First test if (w, father) is in the matching, if get back in the tree, until we find two successive
 * non matching edges, or until we arrive at u.
 *
 * return 1 if (u,v) is valid, 0 otherwise.
 */
int spasm_valid_edge(spasm *A, spasm *TA, int u, spasm_rc *At, spasm_rc *mark, spasm_rc *match, int r){
  int *Ap, *tAp, *Aj, *tAj, *Atr, *Atc, *markr, *markc, *matchr, *matchc;
  int w, px, father, ok;
  
  assert(A != NULL);
  assert(TA != NULL);
  assert(At != NULL);

  Ap = A->p;
  tAp = TA->p;
  Aj = A->j;
  tAj = TA->j;
  Atr = At->r;
  Atc = At->c;
  markr = mark->r; // rows that are upper than u in the tree.
  markc = mark->c; // cols that are upper than u in the tree.
  matchr = match->r;
  matchc = match->c;
  

  //test if u has backward edge.
  if(r == 1){ // u is a row index.
    for(px = Atr[u]; px < Ap[u+1]; px++){
      w = Aj[px]; // get all critical nodes.

      if(markc[w] != -1){ // w has already been considered.
	continue;
      }

      father = tAj[tAp[w]]; // the father of w is the first entry of the col (w isn't root).

      if(matchc[w] != father){ // if c isn't match with its father, then there is no alterning cycle.
	continue;
      }
      
      // go up in the tree.
      ok = spasm_go_up_in_tree(A, TA, father, u, match, r);
      markc[matchr[u]] = ok;
      markr[u] = ok; //if u, and matchr[u] are valid: 1, 0 otherwise.
      if(ok != 1){ // alterning cycle has been found.
	printf("Alterning cycle edge: (%d, %d)\n", u, matchr[u]);
	matchc[matchr[u]] = -1; // get u and v out of the matching.
	matchr[u] = -1;
	return 0;
      }
      	
    }

      
  }else{ //u is a column.
    for(px = Atc[u]; px < tAp[u+1]; px++){
      w  = tAj[px];
      
      if(markr[w] == 1){ // w is upper than u in the tree, so no alterning cycle.
	continue;
      }

      father = Aj[Ap[w]]; // the father of w is the first entry of the row (w isn't root).

      if(matchc[w] != father){ // if c isn't match with its father, then there is no alterning cycle.
	continue;
      }
      
      // go up in the tree.
      ok = spasm_go_up_in_tree(A, TA, father, u, match, r);
      markr[matchc[u]] = ok;
      markc[u] = ok;
      if(ok != 1){ // alterning cycle has been found.
	printf("Alterning cycle edge: (%d, %d)\n", u, matchc[u]);
	matchr[matchc[u]] = -1; // get u and v out of the matching.
	matchc[u] = -1;
	return 0;
      }
      
    }
  }

  return 1;
  
}


/*
 *perform the DFS to create the tree and the matching. Then, for each edge of the matching,
 *check if there is a backward egde on the tree. If so, test if there is an alterning cycle.
 *If so, remove the closest edge from the root from the matching.
 */
spasm_rc *spasm_ur_matching(spasm *A){
  spasm_rc *match;
  spasm *TA;
  spasm_rc *At, *first_passage, *mark;
  int i, j, count, count_col =0;
  int *Atr, *Atc, *fpr, *fpc, *matchr, *matchc, *markr, *markc, *col_TO;

  //allocate memory.
  match = spasm_rc_alloc(A->n, A->m);

  //get workspace.
  first_passage = spasm_rc_alloc(A->n, A->m);
  At = spasm_rc_alloc(A->n, A->m);
  mark = spasm_rc_alloc(A->n, A->m);
  TA = spasm_transpose(A, 0);
  col_TO = spasm_malloc(A->m * sizeof(int));

  Atr = At->r;
  Atc = At->c;
  fpr = first_passage->r;
  fpc = first_passage->c;
  markr = mark->r;
  markc = mark->c;
  matchr = match->r;
  matchc = match->c;

  //init workspace.
  for(i = 0; i < A->n; i++){
    Atr[i] = A->p[i]; // no node in the tree yet.
    fpr[i] = -1; // no row seen yet.
    matchr[i] = -1; // no edge in the matching yet.
    markr[i] = -1;
  }

  for(i = 0; i < A->m; i++){
    Atc[i] = TA->p[i];
    fpc[i] = -1;
    matchc[i] = -1;
    markc[i] = -1;
  }

  //DFS:
  count = spasm_dfs_bipartite(A, TA, 0, first_passage, match, At, col_TO, 0, &count_col);

  while(count != (A->n + A->m)){ // not a tree but a forest.
    for(i = 0; i < A->n; i++){
      
      if(fpr[i] == -1){ // not marked yet.
	break;
      }
    }

    if(i == A->n){ // Empty cols not seen by process
      break;
    }
    
    count = spasm_dfs_bipartite(A, TA, i, first_passage, match, At, col_TO, count, &count_col);
   }

  //count matching edge.
  count = 0;
  for(i = 0; i < count_col; i++){
    if(matchc[col_TO[i]] != -1){
      count++;
    }
  }

  fprintf(stderr,"nb edges in matching: %d\n", count);


  /* for(i = 0; i <A->m; i++){ */
  /*   printf("col_TO[%d] = %d\n", i, col_TO[i]); */
  /* } */

  //go up phase:
  count = 0;
  for(i = 0; i < count_col; i++){ // look at the column. no root.
    
    if(matchc[col_TO[i]] == -1){ // col i is not in the matching.
      continue;
    }

    j = col_TO[i];
    if(matchc[j] == TA->p[j]){ // the first row on col j is the father of j (no root).
      //j is match with its father.

      // test whether the father of i is in an alterning cycle.
      count += spasm_valid_edge(A, TA, TA->p[j], At, mark, match, 1);
     
    }else{ //i is matched with one of its sons
      count += spasm_valid_edge(A, TA, j, At, mark, match, 0);
      
    }

    

  }
  

  fprintf(stderr,"number of valid edgde: %d\n", count);
  
  

  //free workspace.
  spasm_rc_free(first_passage);
  spasm_rc_free(At);
  spasm_rc_free(mark);
  spasm_csr_free(TA);
  free(col_TO);

  //return.
  return match;
}
