#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(){

  spasm *A, *U;
  super_spasm *L;
  spasm_triplet *T1, *T2;
  spasm_lu *LU;
  int r, n, m, i, j, unz, lnz, prime, li, ui, top, npiv, found;
  int *xi, *yi, *qinv, *R, *Lp, *Up, *w;
  spasm_GFp *x, *y;


  /*loading matrix */
  T1 = spasm_load_sms(stdin, 42013);
  T2 = spasm_triplet_alloc(T1->m, T1->n, T1->nzmax, T1->prime, 1);
  for(i = 0; i < T1->nz; i++){
    if(T1->i[i] & 0x1){
      T2->i[i] = T1->i[i];
      T2->j[i] = T2->j[i];
      T2->x[i] = T2->x[i];
    }
  }
  spasm_triplet_free(T1);
  A = spasm_compress(T2);
  spasm_triplet_free(T2);

  n = A->n;
  m = A->m;

  r = spasm_min(n, m);
  prime = A->prime;

  /* Allouer U et L */
  unz = 4 * spasm_nnz(A) + m;
  lnz = 4 * spasm_nnz(A) + m; // educated gess.

  U = spasm_csr_alloc(r, m, unz, prime, 1);
  L = super_spasm_alloc(n, n, n, lnz, prime, 1);

  Up = U->p;
  Lp = L->M->p;

  /* Get workspace */
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(xi, 3*m);
  spasm_vector_zero(x, m);
  qinv = spasm_malloc(m * sizeof(int));
  R = spasm_malloc(r * sizeof(int));
  // w = spasm_malloc(n * sizeof(int));
  // spasm_vector_zero(w, n);

  /* initialize workspace */
  for(i = 0; i < m; i++){
    qinv[i] = -1;
  }

  for(i = 0; i < r; i++){
    R[i] = -1;
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
    if (unz + m > U->nzmax) {
      spasm_csr_realloc(U, 2 * U->nzmax + m);
    }

    /* triangular solve */
    top = spasm_sparse_forward_solve(U, A, i, xi, x, qinv);

    /* search for pivot */
    found = super_spasm_find_pivot(xi, x, top, U, L, &unz, &lnz, li, ui, i, qinv, R);

    li++;
    ui += found;
    npiv += found;
    // w[i] = found; // w[i] : nombre de pivots sur la ligne i.

  }

/* Finalize L and U */

  Up[ui] = unz;
  spasm_csr_resize(U, ui, m);
  spasm_csr_realloc(U, -1);

  Lp[li] = lnz;
  spasm_csr_resize(L->M, li, n);
  spasm_csr_realloc(L->M, -1);

  //free(x);
  //free(xi);

  /* Check result */
  LU = spasm_LU(A, NULL, 1);

  assert(npiv == LU->U->n);
  r = npiv;

  // check qinv :
  for(i = 0; i < m; i++){
    if(qinv[i] != LU->qinv[i]){
      printf("Not ok qinv column %d \n", i);
      exit(0);
    }
  }

  // check U :
for (i = 0; i < npiv; i++){
    if(U->p[i+1] != LU->U->p[i+1]){
      printf("Not ok U row %d\n", i);
      exit(0);
    }
    for(j = U->p[i]; j < U->p[i+1]; j++){
      if(U->j[j] != LU->U->j[j]){
	printf("Not ok U row %d colums \n", i);
	exit(0);
      }
      if(U->x[j] != LU->U->x[j]){
	printf("Not ok U row %d entries \n", i);
	exit(0);
	}
    }
  }


/* check L */


 printf("ok super_find_pivot \n");

  /* free memory */
  spasm_csr_free(U);
  spasm_csr_free(A);
  super_spasm_free(L);
  spasm_free_LU(LU);
  free(x);
  free(xi);
  free(qinv);
  free(R);
  //  free(w);

}
