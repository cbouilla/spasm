#include <stdio.h>
#include <assert.h>
#include "spasm.h"

#ifdef SPASM_TIMING
extern int64 reach, scatter, data_shuffling;
#endif


int main() {
  spasm_triplet *T;
  spasm *A, *U, *L;
  spasm_lu *LU;
  int *p;
  double start_time, end_time;

  T = spasm_load_sms(stdin, 42013);
  printf("A : %d x %d with %d nnz (loaded modulo 42013)\n", T->n, T->m, T->nz);
  if (T->n < T->m) {
    printf("[rank] transposing matrix\n");
    spasm_triplet_transpose(T);
  }
  A = spasm_compress(T);
  spasm_triplet_free(T);

  start_time = spasm_wtime();

#ifdef SPASM_SORT_ROWS
  p = spasm_row_sort(A);
#else
  p = spasm_cheap_pivots(A);
#endif

  LU = spasm_LU(A, p, 1);
  end_time = spasm_wtime();
  printf("\n");

  U = LU->U;
  L = LU->L;

  printf("LU factorisation (+ sort took) %.2f s\n", end_time - start_time);
  printf("U :  %d x %d with %d nnz\n", U->n, U->m, spasm_nnz(U));
  printf("L :  %d x %d with %d nnz\n", L->n, L->m, spasm_nnz(L));

#ifdef SPASM_TIMING
  printf("----------------------------------------\n");
  printf("reach   : %12" PRId64 "\n", reach);
  printf("scatter : %12" PRId64 "\n", scatter);
  printf("misc    : %12" PRId64 "\n", data_shuffling);
#endif
  printf("----------------------------------------\n");
  printf("rank of A = %d\n", U->n);
  spasm_free_LU(LU);


  spasm_csr_free(A);
  return 0;
}
