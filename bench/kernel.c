#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A, *A_t, *K, *A_clean;
  int ch, output, prime, i, n, m, r, h, *row_permutation, *column_permutation;
  double A_density, K_density;

  /* options descriptor */
  struct option longopts[5] = {
    { "verbose",       no_argument,       NULL,    'v' },
    { "tabulated",     no_argument,       NULL,    't' },
    { "matrix",        no_argument,       NULL,    'm' },
    { "modulus",       required_argument, NULL,    'p' },
    { NULL,            0,                 NULL,     0  }
  };

  output = 0;
  prime = 42013;

  while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    switch (ch) {
    case 'p':
      prime = atoi(optarg);
      break;
    case 'v':
      output = 0;
      break;
    case 't':
      output = 1;
      break;
    case 'm':
      output = 2;
      break;
    default:
      printf("Unknown option\n");
      exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  spasm_triplet_free(T);
  n = A->n;
  m = A->m;

  /* remove zero rows */
  int n_cheap;
  row_permutation = spasm_cheap_pivots(A, & n_cheap); /* this stacks zero row at the bottom */
  A_clean = spasm_permute(A, row_permutation, SPASM_IDENTITY_PERMUTATION, SPASM_WITH_NUMERICAL_VALUES);
  spasm_csr_free(A);

   A = A_clean;
  for(i = 0; i < n; i++) {
    if (spasm_row_weight(A, i) == 0) {
      fprintf(stderr, "[kernel] ignoring %d empty rows\n", n - i);
      A->n = i;
      n = i;
      break;
    }
  }
  free(row_permutation);

  A_t = spasm_transpose(A, SPASM_WITH_NUMERICAL_VALUES);
  column_permutation = spasm_cheap_pivots(A_t, &n_cheap);

  K = spasm_kernel(A_t, column_permutation);

  r = K->n;
  K_density = 100.0 * spasm_nnz(K) / (K->n * K->m);
  A_density = 100.0 * spasm_nnz(A) / (A->n * A->m);

  /* identify heaviest row */
  h = -1;
  for (i = 0; i < r; i++) {
    h = spasm_max(h, spasm_row_weight(K, i));
  }

  switch (output) {
  case 0:
    /* output some stats in human-readable form */
    printf("kernel of dimension %d\n", K->n);
    printf("kernel NNZ = %d   /   density = %.2f %%  /   increase wrt A : %.1f\n", spasm_nnz(K), K_density, K_density / A_density);
    printf("heaviest row NNZ = %d   / density = %.2f %%\n", h, 100.0 * h / n);
    break;
  case 1:
    /* output some stats in machine-readable form */
    printf("%d \t %d \t %d\n", K->n,  spasm_nnz(K), h);
    break;

  case 2:
    /* output the kernel basis matrix */
    spasm_save_csr(stdout, K);
    break;
  }

  free(column_permutation);
  spasm_csr_free(A);
  spasm_csr_free(A_t);
  spasm_csr_free(K);
  return 0;
}
