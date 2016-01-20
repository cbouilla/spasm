#include<assert.h>
#include"spasm.h"

/*
 * replace u by u + v.
 * return value is number of non zero entries in new vector u.
 */
int spasm_add_vectors(spasm_GFp *u, int *ui, int unz, spasm_GFp *v, int *vi, int vnz, int size){
  int *w;
  int i, j;

  // get workspace.
  w = spasm_malloc(size * sizeof(int));

  //initialize tables 
  spasm_vector_zero(w, size);

  for(j = 0; j < unz; j++){
    i = ui[j];
    w[i] = 1;
  }

  for(j = 0; j < vnz; j++){
    i = vi[j];
    if(w[i] == 0){ // u[i] = 0
      ui[unz] = i;
      u[i] = v[i]; // <-- 0 + v[i].
      unz++;
    }
    else { //u[i] != 0 , i already in ui.
      u[i] += v[i]; // <-- u[i] + v[i].
    }
  }

  return unz;
}

/*
 * Given a system L, allocate the right-hand side member.
 */
void spasm_system_right_hand_init(spasm_system *L, int nnz){
  spasm *M, *B;
  int Mn, size, prime;
  int *Bp;
  
  //check inputs.
  assert(!L->B);

  M = L->M;
  Mn = M->n; // number of rows of M.
  size = L->rect; // size of the rectangular part of M.
  size += Mn; // size of the right-hand side member
  prime = M->prime;

  // Memorie allocation
  B = spasm_csr_alloc(1, size, nnz, prime, SPASM_WITH_NUMERICAL_VALUES);

  //initialisation
  Bp = B->p;
  Bp[0] = 0;
  Bp[1] = nnz;

  L->B = B;
}

void spasm_first_system_init(spasm_system *L, int i){
  int *Bj;
  spasm_GFp *Bx;

  //Initialise B
  spasm_system_right_hand_init(L, 1);

  Bj = L->B->j;
  Bx = L->B->x;

  Bj[0] = i;
  Bx[0] = 1;

 
}

/*
 * Given the solution x and it pattern xi, dispatchthe coefficient of
 * x for next step.
 */

void spasm_lazy_right_hand_update(int d, int k, spasm_system *L0, spasm_system *L1, int bound, int *xi, int *x, int top, int stop, int **p){
  spasm *B0, *B1;
  int i, i_new, j, d1, d0, start, nz0, nz1;
  int *B0j, *B1j;
  spasm_GFp *B0x, *B1x;

 //check inputs.
  assert(d > 0);
  assert(k >= 0);

  d1 = L0->diag; // diagonale du prochain système sur l'intervalle de lignes k
  d0 = L1->diag; // diagonale du prochain système sur l'intervalle de colonnes k+d

  assert(d1 < d);
  assert(d0 < d);

  start = L0->M->m;

  nz0 = L0->Bnz;
  nz1 = L1->Bnz;

  if(!L0->B){
    int nzmax;
    nzmax = L0->M->n;
    nzmax += L0->rect;
    spasm_system_right_hand_init(L0, nzmax);
  }

  if(!L1->B){
    int nzmax;
    nzmax = L1->M->n;
    nzmax += L1->rect;
    spasm_system_right_hand_init(L1, nzmax);
  }

  B0 = L0->B;
  B1 = L1->B;
  B0j = B0->j;
  B0x = B0->x;
  B1j = B1->j;
  B1x = B1->x;

  if(d1 > 0){
    if(d0 > 0){ // Cas ou les deux systèmes suivants sont sur une diagonale supérieure.
      for(j = top; j < stop; j++){
	i = xi[j];
	if(i < bound){ //On remplit la partie gauche du second membre du système 1.
	  B1j[(nz1)] = i;
	  B1x[(nz1)] = x[i];
	  (nz1)++;
	}
	else { //On remplit la partie droite du second membre du système 0.
	  i_new = i-bound;
	  i_new += start;
	  B0j[(nz0)] = i_new;
	  B0x[(nz0)] = x[i];
	  (nz0)++;
	}
      }
    }
    else { // Cas où le système 0 est sur la diagonale princiaple.
      for(j = top; j < stop; j++){
	i = xi[j];
	if(i < bound){ //On remplit la partie gauche du second membre du système 1.
	  B1j[(nz1)] = i;
	  B1x[(nz1)] = x[i];
	  (nz1)++;
	}
	else { //On remplit la partie droite du second membre du système 0. 
	  i_new = i - bound;
	  i_new += start;
	  i_new = p[k][i_new];
	  B0j[(nz0)] = i_new;
	  B0x[(nz0)] = x[i];
	  (nz0)++;
	}
      }
    }
  }
  else { // Cas où le système 1 est sur la diagonale principale.
    if(d0 > 0){ // Cas où le système 0 est sur une diagonale supérieure.
      for(j = top; j < stop; j++){
	i = xi[j];
	if(i < bound){
	  i_new = p[k+1][i];
	  B1j[(nz1)] = i_new;
	  B1x[(nz1)] = x[i];
	  (nz1)++;
	}
	else { 
	  i_new = i-bound;
	  i_new += start;
	  B0j[(nz0)] = i_new;
	  B0x[(nz0)] = x[i];
	  (nz0)++;
	}
      }   
    }
    else{ // Cas où le système 0 est sur une diagonale principale.
      for(j = top; j < stop; j++){
	i = xi[j];
	if(i < bound){
	  i_new = p[k+1][i];
	  B1j[(nz1)] = i_new;
	  B1x[(nz1)] = x[i];
	  (nz1)++;
	}
	else {
	  i_new = i - bound;
	  i_new += start;
	  i_new = p[k][i_new];
	  B0j[(nz0)] = i_new;
	  B0x[(nz0)] = x[i];
	  (nz0)++;
	}
      }
    }
  }

}

/*
 * given d and k, index, find the good permutation for the system y M_{d,k} = x.
 * bound is the number of right part column
 * pk is a permutation over the row of the k-th row interval.
 * p_new is the new permutation returned by the function
 */
void spasm_new_lazy_permutation(int bound, const int *Lperm, int *p_new, int vec_size){
  int i, diff;

  for(i = 0; i < bound; i++){
    p_new[i] = i;
  }

  diff = vec_size - bound;

  for(i = 0; i < diff; i++){
    p_new[bound+i] = Lperm[i]+bound;
  }
}

/*
 * Given a system L, find out which is the next system corresponding to the "left part" of the solution.
 * return the index of the row interval. 
 */
int spasm_next_left_system(spasm_system **L, int k, int d){
  spasm_system *tmp;
 
  k = k+1;
  if(d == 0){
    return k;
  }

  tmp = L[k];
  while(tmp->diag > d){
    tmp = tmp->next;
  }
  if (tmp->diag == d){
    return k;
  }
  else{
    return spasm_next_left_system(L, k, d-1);
  }
}

/*
 * Solve the system L, and dispatch coefficients for the next step.
 */
void spasm_lazy_system(spasm_system **L, int k, int l, int d, int **p){
  spasm *B, *M;
  int bound, size, top, left, d_new;
  int *xi;
  spasm_GFp *x, *Lp;
  spasm_system *Lk, *Lleft;

  Lk = L[l];

  //check entries
  assert(Lk->B);
  assert(Lk->diag == d);

  bound = Lk->rect;
  left = Lk->left; // next system corresponding to the left part of the solution.

  assert(left > k); // next system is "under" the curent one.

  B = Lk->B;
  size = B->m;
  M = Lk->M;
  Lp = Lk->p;

  // Allocate workspace
  x = spasm_malloc(size * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * size * sizeof(int));

  //initialize vector
  spasm_vector_zero(x, size);
  spasm_vector_zero(xi, size);

  // Solve system
  top = spasm_sparse_backward_solve(M, B, 0, xi, x, Lp, bound);

  // Free right-hand side member.
  free(B);

  // Find diagonal number of next "left" system.
  d_new = d - left;

  assert(d_new < d); // "left system" is on a previous diagonal
  assert(d_new >= 0);

  Lleft = L[l+left];
  while(Lleft->diag > d_new){
    Lleft = Lleft->next;
  }

  // check diagonal number.
  assert(Lleft->diag == d_new);

  // find next "right system" :
  Lk = Lk->next;

  // dispatch coefficient for next systems.
  spasm_lazy_right_hand_update(d, k, Lk, Lleft, bound, xi, x, top, size, p);

  // free x, xi.
  free(x);
  free(xi);

}

/*
 * Given a spasm_list L, a matrix Lnew, an int diag and a permutation p, add Lnew diag and p at the begining of L.
 */
spasm_system * spasm_system_update(spasm_system *L, spasm *M, int *p, int rect, int left, int diag){
  spasm_system *LL;

  LL = spasm_malloc(sizeof(spasm_system));
  LL->M = M;
  LL->B = NULL; // No right-hand side member stocked
  LL->diag = diag;
  LL->rect = rect;
  LL->left = left;
  LL->p = p;
  LL->next = L;

  return LL;
}

/*
 * Free a spasm_system L.
 */
spasm_system * spasm_system_clear(spasm_system *L){

  if(L == NULL){
    return NULL;
  }
  else {
    spasm_system *tmp;
    tmp = L->next;
    spasm_csr_free(L->M);
    if(L->B != NULL){
      spasm_csr_free(L->B);
    }
    free(L->p);
    free(L);
    return spasm_system_clear(tmp);
  }
}


/*
 * Given a diagonal d and index k and i, do the lazy computation:
 * solving successive systems, and final product.
 * The return value is ynz, the number of non-zero entries in the
 * output vector.
 */
int spasm_lazy_computation(int d, int k, int i, spasm_system **S, spasm_GFp *u, int *ui, int usize, spasm **A, int **p){
  int diag, l, nsys, nnz, nztmp;
  spasm_system **L;
  int *xi;
  spasm_GFp *x;
  // spasm *B;

  //check inputs
  assert(d > 1);
  if(A == NULL){
    return 0;
  }

  //Get workspace.
  L = spasm_malloc(d * sizeof(spasm_system*));

  // Initialise L.
  for(l = 0; l < d; l++){
    L[l] = S[k+l];

    assert(L[l]->diag <= d); // check systems.
  }

  // Initialise variable.
  diag = d; // Current diagonal
  nsys = 1;

  //Initialise the first system.

  spasm_first_system_init(L[0], i);

  while(diag > 1){ 
    diag = diag - 1; //previous diagonal

    for(l = 0; l < nsys; l++){
      if(L[l]->diag == diag){ // if L_{k+l, diag} exist.

	//solve the system L and dispatch solution coefficient for next step.
	spasm_lazy_system(L, k, l, diag, p);

	// update L for next step.
	L[l] = L[l]->next;
      }
    }
    for(l = nsys; l < d; l++){
      if(L[l]->diag == diag) L[l] = L[l]->next; // Update L for next step.
    }
    nsys++;
  }

  // Last diagonal system solving and vector-matrix product.
  
  // Get workspace.
  x = spasm_malloc(usize * sizeof(spasm_GFp));
  xi = spasm_malloc(usize * sizeof(int));

  //Initialize x, xi.
  spasm_vector_zero(x, usize);
  spasm_vector_zero(xi, usize);

  nnz = 0;

  for(l = 0; l < nsys; l++){

    assert(L[l]->diag == 0); //check diagonal number.
    // B = L[l]->B;

    if(L[l]->B){
      nztmp = spasm_solve_and_product(L[l]->M, A[l], L[l]->B, 0, x, xi, p[k+l]); //last system solving + product.
      nnz = spasm_add_vectors(u, ui, nnz, x, xi, nztmp, usize);
      free(L[l]->B);
      
    }

  }

  //free workspace
  free(x);
  free(xi);

  //return number of non-zero entries in u.

  return nnz;

}

