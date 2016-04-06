#include <stdio.h>
#include <assert.h>
#include <getopt.h>

#include "spasm.h"

#ifdef SPASM_TIMING
extern int64 reach, scatter, data_shuffling;
#endif
#include <signal.h>
#include <setjmp.h>


jmp_buf Env;

void alarm_handler(int dummy) {
  longjmp(Env, 1);
}



int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *U, *L;
  spasm_lu *LU;
  int r, n, m, *p, ch, prime, allow_transpose, sort_strategy, keep_L, timer;
  double start_time, end_time;

  prime = 42013;
  sort_strategy = 1; // cheap pivots by default
  allow_transpose = 1;
  keep_L = 0;
  timer = -1;

  /* options descriptor */
  struct option longopts[7] = {
    { "sort-rows",    no_argument,       NULL,          's' },
    { "keep-rows",    no_argument,       NULL,          'k' },
    { "no-transpose", no_argument,       NULL,          'a' },
    { "modulus",      required_argument, NULL,          'p' },
    { "keep-L" ,      no_argument,       NULL,          'l' },
    { "max-time",     required_argument, NULL,          't' },
    { NULL,           0,                 NULL,           0  }
  };

  while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    switch (ch) {
    case 'p':
      prime = atoi(optarg);
      break;
    case 'k':
      sort_strategy = 0;
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
    case 't':
      timer = atoi(optarg);
      break;
    default:
      printf("Unknown option\n");
      exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  T = spasm_load_sms(stdin, prime);
  fprintf(stderr, "A : %d x %d with %d nnz (density = %.3f %%) -- loaded modulo %d\n", T->n, T->m, T->nz, 100.0 * T->nz / (1.0 * T->n * T->m), prime);
  if (allow_transpose && (T->n < T->m)) {
    fprintf(stderr, "[rank] transposing matrix : ");
    fflush(stderr);
    start_time = spasm_wtime();
    spasm_triplet_transpose(T);
    fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
  }
  A = spasm_compress(T);
  spasm_triplet_free(T);

  start_time = spasm_wtime();

  switch(sort_strategy) {
  case 0:
    p = NULL;
    break;
  case 1:
    fprintf(stderr, "[rank] finding cheap pivots : ");
    fflush(stderr);
    start_time = spasm_wtime();
    p = spasm_cheap_pivots(A);
    fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
    break;
  case 2:
    fprintf(stderr, "[rank] sorting rows : ");
    fflush(stderr);
    start_time = spasm_wtime();
    p = spasm_row_sort(A);
    fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
    break;
  }

  if (timer > 0) {
    signal(SIGALRM, alarm_handler);
    alarm(timer);

    if (setjmp(Env) != 0) {
      fprintf(stderr, "\nTimeout after %d seconds\n", timer);
      exit(2);
    }
  }

  LU = spasm_LU(A, p, keep_L);
  end_time = spasm_wtime();

  U = LU->U;
  r = U->n;
  n = A->n;
  m = A->m;

  printf("%.2f\n", end_time - start_time);
  fprintf(stderr, "U :  %d x %d with %d nnz (density = %.1f %%)\n", r, m, spasm_nnz(U), 100.0 * spasm_nnz(U) / (1.0*r*m - r*r/2.0));
  if (LU->L != NULL) {
    L = LU->L;
    fprintf(stderr, "L :  %d x %d with %d nnz (density =%.1f %%)\n", L->n, r, spasm_nnz(L), 100.0 * spasm_nnz(L) / (1.0*r*n - r*r/2.0));
  }

#ifdef SPASM_TIMING
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "reach   : %12" PRId64 "\n", reach);
  fprintf(stderr, "scatter : %12" PRId64 "\n", scatter);
  fprintf(stderr, "misc    : %12" PRId64 "\n", data_shuffling);
  fprintf(stderr, "----------------------------------------\n");
#endif
  // printf("rank of A = %d\n", U->n);
  spasm_free_LU(LU);


  spasm_csr_free(A);
  return 0;
}
