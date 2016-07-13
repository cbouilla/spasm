/* indent -nfbs -i2 -nip -npsl -di0 -nut rank_gplu.c */
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
  dummy++;
  longjmp(Env, 1);
}



int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *U, *L;
  spasm_lu *LU;
  int r, n, m, *p, ch, prime, allow_transpose, sort_strategy, keep_L, timer;
  double start_time, end_time;

  prime = 42013;
  sort_strategy = 1;            /* cheap pivots ON by default */
  allow_transpose = 1;          /* transpose ON by default */
  keep_L = 0;
  timer = -1;

  /* options descriptor */
  struct option longopts[7] = {
    {"sort-rows", no_argument, NULL, 's'},
    {"keep-rows", no_argument, NULL, 'k'},
    {"no-transpose", no_argument, NULL, 'a'},
    {"modulus", required_argument, NULL, 'p'},
    {"keep-L", no_argument, NULL, 'l'},
    {"max-time", required_argument, NULL, 't'},
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
  if (allow_transpose && (T->n < T->m)) {
    fprintf(stderr, "[rank] transposing matrix : ");
    fflush(stderr);
    start_time = spasm_wtime();
    spasm_triplet_transpose(T);
    fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
  }
  A = spasm_compress(T);
  spasm_triplet_free(T);
  n = A->n;
  m = A->m;


  start_time = spasm_wtime();

  switch (sort_strategy) {
  case 0:
    p = NULL;
    break;

  case 1:
    fprintf(stderr, "[rank] finding cheap pivots : ");
    fflush(stderr);
    int n_cheap, i, j, k;
    int *Ap = A->p;
    int *Aj = A->j;
    int *qinv = spasm_malloc(m * sizeof(int));
    int *p = spasm_malloc(n * sizeof(int));
    n_cheap = spasm_find_pivots(A, p, qinv);

    /* build qinv to reflect the changes in p */
    for (j = 0; j < m; j++) {
      qinv[j] = -1;
    }

    /* pivotal column first, in row-order */
    k = 0;
    for (i = 0; i < n_cheap; i++) {
      j = Aj[Ap[p[i]]];         /* the pivot is the first entry of each row */
      qinv[j] = k++;
    }

    /* put remaining non-pivotal columns afterwards, in any order */
    for (j = 0; j < m; j++) {
      if (qinv[j] == -1) {
        qinv[j] = k++;
      }
    }

    spasm *B = spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
    free(p);
    free(qinv);
    p = NULL;

    spasm_csr_free(A);
    A = B;

    fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
    break;

  case 2:
    fprintf(stderr, "[rank] sorting rows : ");
    fflush(stderr);
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

  printf("%.2f\n", end_time - start_time);
  fprintf(stderr, "U :  %d x %d with %d nnz (density = %.1f %%)\n", r, m, spasm_nnz(U), 100.0 * spasm_nnz(U) / (1.0 * r * m - r * r / 2.0));
  if (LU->L != NULL) {
    L = LU->L;
    fprintf(stderr, "L :  %d x %d with %d nnz (density =%.1f %%)\n", L->n, r, spasm_nnz(L), 100.0 * spasm_nnz(L) / (1.0 * r * n - r * r / 2.0));
  }
#ifdef SPASM_TIMING
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, "reach   : %12" PRId64 "\n", reach);
  fprintf(stderr, "scatter : %12" PRId64 "\n", scatter);
  fprintf(stderr, "misc    : %12" PRId64 "\n", data_shuffling);
  fprintf(stderr, "----------------------------------------\n");
#endif

  printf("rank of A = %d\n", U->n);

  free(p);
  spasm_free_LU(LU);
  spasm_csr_free(A);
  return 0;
}
