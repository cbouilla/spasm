#include <stdio.h>
#include <assert.h>
#include "spasm.h"

#ifdef SPASM_TIMING
extern int64 reach, scatter, data_shuffling;
#endif


int main() {
  spasm_triplet *T;
  spasm *A, *L, *U;
  spasm_lu *LU;
  double start_time, end_time;
  int *p;

  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  printf("A : %d x %d with %d nnz\n", A->n, A->m, spasm_nnz(A));

  start_time = spasm_wtime();

#ifdef SPASM_SORT_ROWS
  p = spasm_row_sort(A);
#else
  p = SPASM_IDENTITY_PERMUTATION;
#endif

  LU = spasm_LU(A, p, 0);
  end_time = spasm_wtime();
  printf("\n");

  //  L = LU->L;
  U = LU->U;

  printf("LU factorisation (+ sort took) %.2f s\n", end_time - start_time);
  //  printf("L :  %d x %d with %d nnz\n", L->n, L->m, spasm_nnz(L));
  printf("U :  %d x %d with %d nnz\n", U->n, U->m, spasm_nnz(U));
  spasm_free_LU(LU);

#ifdef SPASM_TIMING
  printf("reach   : %12" PRId64 "\n", reach);
  printf("scatter : %12" PRId64 "\n", scatter);
  printf("misc    : %12" PRId64 "\n", data_shuffling);
#endif

  spasm_csr_free(A);
  return 0;
}
