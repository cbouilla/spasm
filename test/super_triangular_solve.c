#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(){
  spasm *A;
  spasm_triplet *T1, *T2;
  super_spasm *L, *U;
  int r, n, m, i, j, unz, lnz, prime, li, ui, top, npiv, found;
  int *xi, *yi, *qinv, *Lp, *Up;
  spasm_GFp *x, *y, *u;

  /* ---------- creat matrix L --------------*/
  /* Given a matrix A, we compute it's LU
   * decomposition using super_spasm_find_pivot
   * function, in order to get L super_spasm.
   */

 /*loading matrix A */
  T1 = spasm_load_sms(stdin, 42013);
  T2 = spasm_triplet_alloc(T1->n, T1->m, T1->nzmax, T1->prime, 1);
  j = 0;
  for(i = 0; i < T1->nz; i++){
    if(T1->i[i] & 0x1){
      T2->i[j] = T1->i[i];
      T2->j[j] = T1->j[i];
      T2->x[j] = T1->x[i];
      j++;
    }
  }
  T2->nz = j;
  spasm_triplet_free(T1);
  A = spasm_compress(T2);
  spasm_triplet_free(T2);


  /* super sparse decomposition LU of A */
  n = A->n;
  m = A->m;

  r = spasm_min(n, m);
  prime = A->prime;

  /* Allouer U et L */
  unz = 4 * spasm_nnz(A) + m;
  lnz = 4 * spasm_nnz(A) + m; // educated gess.

  U = super_spasm_alloc(n, r, m, unz, prime, 1);
  L = super_spasm_alloc(n, n, n, lnz, prime, 1);

  Up = U->M->p;
  Lp = L->M->p;

  /* Get workspace */
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(xi, 3*m);
  spasm_vector_zero(x, m);
  qinv = spasm_malloc(m * sizeof(int));
  // w = spasm_malloc(n * sizeof(int));
  //spasm_vector_zero(w, n);

  /* initialize workspace */
  for(i = 0; i < m; i++){
    qinv[i] = -1;
  }

  for(i = 0; i < r; i++){
    Up[i] = 0;
  }

  for(i = 0; i < n; i++){
    Lp[i] = 0;
  }

  unz = lnz = 0;
  li = ui = 0;
  npiv = 0;

  /* main loop : compute L[i] and U[i] */
  for(i = 0; i < n; i++){
    if(!(i & 0x1)){
      continue;
    }
    
    Lp[li] = lnz;
    Up[ui] = unz;

 /* not enough room in L/U ? realloc twice the size */
    if (lnz + m > L->M->nzmax) {
      spasm_csr_realloc(L->M, 2 * L->M->nzmax + m);
    }
    if (unz + m > U->M->nzmax) {
      spasm_csr_realloc(U->M, 2 * U->M->nzmax + m);
    }

    /* triangular solve */
    top = spasm_sparse_forward_solve(U->M, A, i, xi, x, qinv);

    /* search for pivot */
    found = super_spasm_find_pivot(xi, x, top, U, L, &unz, &lnz, li, ui, i, qinv);

    li++;
    ui += found;
    npiv += found;
    // w[i] = found; // w[i] : nombre de pivots sur la ligne i.

  }

/* Finalize L and free U */
  super_spasm_free(U);

  Lp[li] = lnz;
  spasm_csr_resize(L->M, li, n);
  spasm_csr_realloc(L->M, -1);
  Lp = L->p;

  free(x);
  free(xi);

  /* -------- Solve a system x*L = y --------- */
  
  /* get workspace */
  x = spasm_malloc(n * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * n * sizeof(int));
  y = spasm_malloc(n * sizeof(spasm_GFp));
  yi = spasm_malloc(n * sizeof(int));
  u = spasm_malloc(n * sizeof(spasm_GFp));

  /* initialize workspace */
  for(i = 0; i < n; i++){
    yi[i] = 0;
    y[i] = 0;
    x[i] = 0; // clear workspace.
    u[i] = 0;
  }

  // rhs pattern :
  yi[0] = 0;
  yi[1] = n/2;
  yi[2] = yi[1] + 1;
  yi[3] = n - 1;

  // rhs values :
  for(i = 0; i < 4; i++){
    y[yi[i]] = i+1;
  }

  top = super_spasm_sparse_solve(L, y, yi, 0, x, xi);

  /* --- check result -------- */
  super_sparse_gaxpy_dense(L, x, u); // u <- x*(L-I)
  for(i = 0; i < n; i++){
    u[i] = (u[i] + x[i]) % prime;
  }

  for(i = 0; i < n ; i++){
    if(u[i] != y[i]){
      //printf("u[%d] = %d : y[%d] = %d\n", i, u[i], i, y[i]);
       printf("not ok col %d \n", i);
    }
  }

  printf("ok super triangular solve \n");

  /* free workspace */
  spasm_csr_free(A);
  super_spasm_free(L);
  free(x);
  free(xi);
  free(y);
  free(yi);
  free(u);
  free(qinv);
}
