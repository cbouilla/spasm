#include <stdio.h>
#include <assert.h>
#include "spasm.h"


int main() {
  spasm_triplet *T;
  spasm *A, *L, *U;
  spasm_lu *LU;
  double start_time, end_time;
  int *p;

  T = spasm_load_triplet(stdin, 7);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  printf("A : %d x %d with %d nnz\n", A->n, A->m, spasm_nnz(A));

  start_time = spasm_wtime();
  LU = spasm_LU(A, SPASM_IDENTITY_PERMUTATION);
  end_time = spasm_wtime();

  L = LU->L;
  U = LU->U;

  printf("LU factorisation took %.2f s\n", end_time - start_time);
  printf("L :  %d x %d with %d nnz\n", L->n, L->m, spasm_nnz(L));
  printf("U :  %d x %d with %d nnz\n", U->n, U->m, spasm_nnz(U));
  spasm_free_LU(LU);


  start_time = spasm_wtime();
  p = spasm_row_sort(A);
  LU = spasm_LU(A, p);
  end_time = spasm_wtime();

  L = LU->L;
  U = LU->U;

  printf("LU factorisation + sort took %.2f s\n", end_time - start_time);
  printf("L :  %d x %d with %d nnz\n", L->n, L->m, spasm_nnz(L));
  printf("U :  %d x %d with %d nnz\n", U->n, U->m, spasm_nnz(U));
  spasm_free_LU(LU);


  spasm_csr_free(A);
  return 0;
}
