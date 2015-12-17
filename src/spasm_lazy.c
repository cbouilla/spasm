#include<assert.h>
#include"spasm.h"

/*
 * Given a diagonal index d>0, index k and i, 
 * permutations p and R, pivots numbers indicators 
 * ri and rj, return 0 if i is in the right side of the vector
 * and 1 otherwise. Update index i for
 * the next stage (solving y_k L_{d-1,k} = x_k)
 * and the position of i in vector .
 */
int spasm_lazy_vector_update(int d, int k, int n_blocks, int *i_ptr, const int **ri, const int **rj, const int **p){
  int i, d_new, dr, bound, l; 

  //check inputs.
  assert(d > 0);
  assert(k >= 0);
  assert(k+1 < n_blocks-d);

  i = *i_ptr;
  bound = rj[k+d][d-1];

  if(d > 1){
    if(i < bound) {
      l = 1;
    }
    else {
      d_new = d-1;
      i = i - bound;
      dr = ri[k][d_new] - ri[k][d_new];
      i = i + rj[k+d][d_new-1];
      i = i + dr;
      *i_ptr = i;
      l = 0;
    }
  }
  else {
    if (i < bound){
      k++;
      i = p[k][i];
      *i_ptr = i;
      l = 1;
    }
    else {
      i = i - bound;
      i = i + ri[k][0];
      i = p[k][i];
      *i_ptr = i;
      l = 0;
    }
  }
  return l;
}


/* 
 * Given rj and ri and a1 for block (k, k+d), find size of vector x_k at stage d.
 */
int spasm_lazy_vector_size(int ri, int rj, int a1){
  int size;
  size = a1 - ri;
  size += rj;

  return size;
}


/*
 * Initialise the right side member given a vector x and its pattern.
 */
spasm *spasm_convert_vector_to_matrix(int *x, int *xi, int size, int prime, int xnz){
  spasm *B;
  int *Bp, *Bj, *Bx, j;

  B = spasm_csr_alloc(size, 1, xnz, prime, 1);
  Bp = B->p;
  Bj = B->j;
  Bx = B->x;

  Bp[1] = xnz;

  for(j = 0; j < xnz; j++){
    Bj[j] = xi[j];
    Bx[j] = x [xi[j]];
  }

  return B;
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
 * Given d a diagonal number and index k, and i, 
 * do all the step of the Lazy computation : 
 * solving equations and final product.
 * The return value is ynz, the number of non-zero
 * entries in the output vector.
 */
int spasm_lazy_product(int d, int i, int k, int N, spasm ***L, const spasm *M, const int **ri, const int **rj, const int **p, const int *a1, const int ***Lperm, spasm_GFp *u, int *ui){
  spasm **B;
  int n_vec, l, j, top, l_new, k_new, Bn, index, Mn, i_new, vnz, unz, prime;
  int *p_new, *count, *vec_size, *vi;
  int **x, **y;
  spasm_GFp **xi, **yi;
  spasm_GFp *v;

  //check inputs
  assert(d > 1);
  assert(k >= 0);
  assert(k+1 < N-d);

  if(L == NULL || M == NULL) return 0;

  Mn = M->n;
  prime = M->prime;

  /* Get workspace */
  B = spasm_malloc(d * sizeof(spasm*));
  xi = spasm_malloc(d * sizeof(int*));
  x = spasm_malloc(d * sizeof(spasm_GFp*));
  yi = spasm_malloc(d * sizeof(int*));
  y = spasm_malloc(d * sizeof(spasm_GFp*));

  v = spasm_malloc(Mn * sizeof(spasm_GFp));
  vi = spasm_malloc(Mn * sizeof(int));
  spasm_vector_zero(v, Mn);
  spasm_vector_zero(vi, Mn);
 

  count = spasm_malloc(d * sizeof(int));
  vec_size = spasm_malloc(d * sizeof(int));
  // find size of the first vector ei

  vec_size[0] = spasm_lazy_vector_size(ri[k][d-1], rj[k+d][d-1], a1[k]);

  assert(i < vec_size[0]);  

  x[0] = spasm_malloc(vec_size[0] * sizeof(spasm_GFp));
  xi[0] = spasm_malloc(1 * sizeof(int));
  y[0] = spasm_malloc(vec_size[0] * sizeof(spasm_GFp));
  yi[0] = spasm_malloc(3 * vec_size[0] * sizeof(int));

  /* Initialise variable and vector */
  spasm_vector_zero(x[0], vec_size[0]);
  spasm_vector_zero(y[0], vec_size[0]);
  spasm_vector_zero(yi[0], 3 * vec_size[0]);
 
  x[0][i] = 1;
  xi[0][0] = i;

  B[0] = spasm_convert_vector_to_matrix(x[0], xi[0], vec_size[0], prime, 1);
  free(x[0]);
  free(xi[0]);
  n_vec = 1;
  vnz = 0;
  
  spasm_vector_zero(count, d);


  /* For each diagonal greater than 1 */
  while(d > 1){
    d = d - 1;

    /* Get next x, xi workspace */
    for(l = 0; l < n_vec +1; l++){
      k_new = k+l;
    vec_size[l] = spasm_lazy_vector_size(ri[k_new][d-1], rj[k_new+d][d -1], a1[k_new]);
    x[l] = spasm_malloc(vec_size[l] * sizeof(spasm_GFp));
    xi[l] = spasm_malloc(vec_size[l] * sizeof(int));
    spasm_vector_zero(xi[l], vec_size[l]);
    spasm_vector_zero(x[l], vec_size[l]);
    }
  /* For each small vector */
    for(l = 0; l < n_vec; l++){
      k_new = k + l;
      Bn = B[l]->n;

  /* Find the good permutation */
      p_new = spasm_malloc(Bn * sizeof(int));
      spasm_new_lazy_permutation(rj[k_new + d][d - 1], Lperm[d-1][k_new], p_new, Bn);

  /* Solve the triangular system */
      top = spasm_sparse_backward_solve(L[d][k_new], B[l], 0, yi[l], y[l], p_new, rj[k_new + d][d -1]);

      //free extra worskpace
      spasm_csr_free(B[l]);
      free(p_new);

  /* For each j in y[l] pattern, dispatch j for next step*/
      for(j = top; j < Bn; j++){
	//find the corresponding vector and the index in it :
	index = yi[l][j];

	l_new = spasm_lazy_vector_update(d, k_new, N, &index, ri, rj, p);
	l_new = l_new + l;

	xi[l_new][count[l_new]] = index;
	count[l_new]++;
	x[l_new][index] = y[l][yi[l][j]];
      }
    }
    
    n_vec++;

    for(l = 0; l < n_vec; l++){
    /* Update matrix B */
      B[l] = spasm_convert_vector_to_matrix(x[l], xi[l], vec_size[l], prime, count[l]); 
      free(x[l]);
      free(xi[l]);

    /*Initialise y, yi for next step */
      k_new = k+l;
      y[l] = spasm_realloc(y[l], vec_size[l] * sizeof(int));
      yi[l] = spasm_realloc(yi[l], 3*vec_size[l] * sizeof(int));
      spasm_vector_zero(y[l], vec_size[l]);
      spasm_vector_zero(y[l], 3*vec_size[l]);

    }

  }

  assert(d == 1);
  /* solving last system and dispatch coeff in v */
  for(l = 0; l < n_vec; l++){
    k_new = k + l;
    Bn = B[l]->n;

    /* Find the good permutation */
    p_new = spasm_malloc(Bn * sizeof(int));
    spasm_new_lazy_permutation(rj[k_new + d][d - 1], Lperm[d-1][k_new], p_new, Bn);

    /* Solve the triangular system */
    top = spasm_sparse_backward_solve(L[d][k_new], B[l], 0, yi[l], y[l], p_new, rj[k_new + d][d - 1]);

    //free extra worskpace
    spasm_csr_free(B[l]);
    free(p_new);

    /* dispatch yi in vi and y in v. */
    for(j = top; j < Bn; j++){
      index = yi[l][j];
      i_new = index + a1[k_new];
      vi[vnz] = i_new;
      v[i_new] = y[l][index];
      vnz ++;
    } 
    //free extra workspace
    free(yi[l]);
    free(y[l]);
  }

  /*calculate product*/
  unz = spasm_sparse_vector_matrix_prod(M, v, vi, vnz, u, ui);
  
  /*free workspace */
  free(v);
  free(vi);

  return unz;
}
