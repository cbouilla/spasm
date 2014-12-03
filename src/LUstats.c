#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *L, *U;
  spasm_lu *LU;
  double start_time, end_time;

  T = spasm_load_triplet(stdin, 7);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  printf("A : %d x %d with %d nnz\n", A->n, A->m, spasm_nnz(A));

  start_time = spasm_wtime();
  LU = spasm_LU(A);
  end_time = spasm_wtime();

  L = LU->L;
  U = LU->U;

  printf("LU factorisation took %.2f s\n", end_time - start_time);
  printf("L :  %d x %d with %d nnz\n", L->n, L->m, spasm_nnz(L));
  printf("U :  %d x %d with %d nnz\n", U->n, U->m, spasm_nnz(U));

  spasm_csr_free(A);
  spasm_free_LU(LU);
  return 0;
}
