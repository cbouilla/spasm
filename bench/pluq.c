/* indent -nfbs -i2 -nip -npsl -di0 -nut pluq.c */
#include <stdio.h>
#include <assert.h>
#include <getopt.h>

#include "spasm.h"

#ifdef SPASM_TIMING
extern int64 reach, scatter, data_shuffling;
#endif

/** computes a PLUQ decomposition. U is always saved in a file named U.sms.
 *  If the keep-L option is provided, then L is also saved in L.sms */

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *U, *L;
  spasm_lu *PLUQ;
  int r, n, m, *p, ch, prime, allow_transpose, sort_strategy, keep_L;
  double start_time, end_time;

  prime = 42013;
  sort_strategy = 1;            /* cheap pivots by default */
  allow_transpose = 1;          /* transpose if more columns than rows by
                                 * default */
  keep_L = 0;

  /* options descriptor */
  struct option longopts[6] = {
    {"sort-rows", no_argument, NULL, 's'},
    {"keep-rows", no_argument, NULL, 'k'},
    {"no-transpose", no_argument, NULL, 'a'},
    {"modulus", required_argument, NULL, 'p'},
    {"keep-L", no_argument, NULL, 'l'},
    {NULL, 0, NULL, 0}
  };

  while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    switch (ch) {
    case 'p':
      prime = atoi(optarg);
      break;
    case 'k':
      sort_strategy = 0;
      break;
    case 't':
      break;
    case 'a':
      allow_transpose = 0;
      break;
    case 's':
      sort_strategy = 2;
      break;
    case 'l':
      keep_L = 1;
      break;
    default:
      printf("Unknown option\n");
      exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  T = spasm_load_sms(stdin, prime);
  printf("A : %d x %d with %d nnz (density = %.3f %%) -- loaded modulo %d\n", T->n, T->m, T->nz, 100.0 * T->nz / (1.0 * T->n * T->m), prime);
  if (allow_transpose && (T->n < T->m)) {
    printf("[rank] transposing matrix : ");
    fflush(stdout);
    start_time = spasm_wtime();
    spasm_triplet_transpose(T);
    printf("%.1f s\n", spasm_wtime() - start_time);
  }
  A = spasm_compress(T);
  spasm_triplet_free(T);

  start_time = spasm_wtime();

  p = NULL;

  switch (sort_strategy) {
  case 0:
    break;
  case 1:
    printf("[rank] finding cheap pivots : ");
    fflush(stdout);
    start_time = spasm_wtime();
    int n_cheap;
    p = spasm_cheap_pivots(A, &n_cheap);
    printf("%.1f s\n", spasm_wtime() - start_time);
    break;
  case 2:
    printf("[rank] sorting rows : ");
    fflush(stdout);
    start_time = spasm_wtime();
    p = spasm_row_sort(A);
    printf("%.1f s\n", spasm_wtime() - start_time);
    break;
  }

  PLUQ = spasm_PLUQ(A, p, keep_L);
  end_time = spasm_wtime();
  printf("\n");

  U = PLUQ->U;
  r = U->n;
  n = A->n;
  m = A->m;

  printf("LU factorisation (+ sort took) %.2f s\n", end_time - start_time);
  printf("U :  %d x %d with %d nnz (density = %.1f %%)\n", r, m, spasm_nnz(U), 100.0 * spasm_nnz(U) / (1.0 * r * m - r * r / 2.0));
  if (PLUQ->L != NULL) {
    L = PLUQ->L;
    printf("L :  %d x %d with %d nnz (density =%.1f %%)\n", L->n, r, spasm_nnz(L), 100.0 * spasm_nnz(L) / (1.0 * r * n - r * r / 2.0));
    FILE *f = fopen("L.sms", "w");
    spasm_save_csr(f, L);
    fclose(f);
  }
#ifdef SPASM_TIMING
  printf("----------------------------------------\n");
  printf("reach   : %12" PRId64 "\n", reach);
  printf("scatter : %12" PRId64 "\n", scatter);
  printf("misc    : %12" PRId64 "\n", data_shuffling);
#endif
  printf("----------------------------------------\n");
  printf("rank of A = %d\n", U->n);

  FILE *f = fopen("U.sms", "w");
  spasm_save_csr(f, U);
  fclose(f);
  spasm_free_LU(PLUQ);

  free(p);
  spasm_csr_free(A);
  return 0;
}
