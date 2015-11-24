#include <stdio.h>
#include <assert.h>
#include "spasm.h"


int main(int argc, char **argv){
  spasm *A, *U, *L;
  spasm_triplet *T;
  spasm_lu *LU;
  spasm_GFp *x;
  int *xi, *qinv, *p, *Lp, *Up;
  int test, i, inew, j, unz, lnz, deff, n, m, r, prime, npiv, top; 

  assert(argc > 1);
  test = atoi(argv[1]);

  /*loading matrix */
  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  r = spasm_min(n, m);
  prime = A->prime;
  deff=0;
  npiv = 0;

  /* Guess of the size of L, U */
  lnz = 4 * spasm_nnz(A) + n;
  unz = 4 * spasm_nnz(A) + n;
  
  /* Get workspace */
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(spasm_GFp));
  spasm_vector_zero(xi, 3*m);
  qinv = spasm_malloc(m * sizeof(int));
  p = spasm_malloc(n * sizeof(int));

  /*allocate L, U */
  L = spasm_csr_alloc(n, r, lnz, prime, 1);
  U = spasm_csr_alloc(r, m, unz, prime, 1);

  Lp = L->p;
  Up = U->p;

  for(i = 0; i < m; i++){
    /* clear workspace */
    x[i] = 0;
    qinv[i] = -1;
  }

  for(i = 0; i < n; i++){
    p[i]=i;
  }
  
  for(i = 0; i < r; i++) {
    Up[i] = 0;
  }
  lnz = unz = 0;

  /* --- Main loop : compute L[i] and U[i] ---------*/
  for(i=0; i<r; i++){
    Lp[i] = lnz;
    Up[i - deff] = unz;

    /* not enough room in L/U ? realloc twice the size */
    if (lnz + m > L->nzmax) {
      spasm_csr_realloc(L, 2 * L->nzmax + m);
    }
    if (unz + m > U->nzmax) {
      spasm_csr_realloc(U, 2 * U->nzmax + m);
    }

    //inew = SPASM_IDENTITY_PERMUTATION;
    top = spasm_sparse_forward_solve(U, A, i, xi, x, qinv);
    npiv += spasm_find_pivot(xi, x, top, U, L, &unz, &lnz, i, &deff, qinv, p, n);
    
  }

  /* Check npiv is the rank of A*/
  if(npiv != i - deff){
    printf("Not ok %d npiv == n - defficiency\n", test);
    exit(0);
  }
  
  /* Finalize L and U */

  Up[i-deff] = unz;
  spasm_csr_resize(U, i -deff, m);
  spasm_csr_realloc(U, -1);

  Lp[n] = lnz;
  spasm_csr_resize(L, n, n -deff);
  spasm_csr_realloc(L, -1);

  free(x);
  free(xi);

  /* check L, U, LU decompisition of A */
  LU = spasm_LU(A, SPASM_IDENTITY_PERMUTATION, 1);
  for(i = 0; i < m; i++){
    if(qinv[i] != LU->qinv[i]){
      printf("Not ok %d qinv column %d \n", test, i);
      exit(0);
    }
  }

  
  for (i = 0; i < n; i++){
    if(p[i] != LU->p[i]){
      printf("Not ok %d p row %d\n", test, i);
      printf("p[i] = %d and should be %d \n", p[i], LU->p[i]);
      exit(0);
    }
    if(L->p[i+1] != LU->L->p[i+1]){
      printf("Not ok %d L row %d\n", test, i);
      exit(0);
    }
    for(j = L->p[i]; j < L->p[i+1]; j++){
      if(L->j[j] != LU->L->j[j]){
	printf("Not ok %d L row %d colums \n", test, i);
	exit(0);
      }
      if(L->x[j] != LU->L->x[j]){
	printf("Not ok %d L row %d entries \n", test, i);
	exit(0);
	}
    }
  }

  
for (i = 0; i < npiv; i++){
    if(U->p[i+1] != LU->U->p[i+1]){
      printf("Not ok %d U row %d\n", test, i);
      exit(0);
    }
    for(j = U->p[i]; j < U->p[i+1]; j++){
      if(U->j[j] != LU->U->j[j]){
	printf("Not ok %d U row %d colums \n", test, i);
	exit(0);
      }
      if(U->x[j] != LU->U->x[j]){
	printf("Not ok %d U row %d entries \n", test, i);
	exit(0);
	}
    }
  }

 printf("Ok %d L*U = A\n", test);
}
