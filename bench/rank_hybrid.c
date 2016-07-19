  /* indent -nfbs -i2 -nip -npsl -di0 -nut rank_hybrid.c */
#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/** Computes the rank of the input matrix using the hybrid strategy */

int main(int argc, char **argv) {

  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, n_times, i, rank, npiv, n, m;
  double start_time, end_time;
  spasm_triplet *T;
  spasm *A, *B;
  int *p, *qinv;
  double density;
  char nnz[6];

  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  /*printf("I am going to FREE\n");
  fflush(stdout);*/
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  spasm_human_format(spasm_nnz(A), nnz);
  fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, nnz);
  
  p = spasm_malloc(n * sizeof(int));
  qinv = spasm_malloc(m * sizeof(int));

  /* 3 iterations of Schur complement, by default */
  n_times = 3;
  if (argc > 1) {
    n_times = atoi(argv[1]);
  }

  start_time = spasm_wtime();
  rank = 0;  
  npiv = spasm_find_pivots(A, p, qinv);
  density = spasm_schur_probe_density(A, p, qinv, npiv, 100);
  
  for (i = 0; i < n_times; i++) {
    int nnz = density * (n - npiv) * (m - npiv);
    char tmp[6];
    spasm_human_format(sizeof(int)*(n-npiv+nnz) + sizeof(spasm_GFp)*nnz, tmp);
    fprintf(stderr, "round %d. Schur complement is %d x %d, estimted density : %.2f (%s byte)\n", i, n-npiv, m-npiv, density, tmp);
  
    if (density > 0.1 || m-npiv < 0.1 * m) {
      break;
    }

    /* compute schur complement, update matrix */
    B = spasm_schur(A, p, npiv);
    spasm_csr_free(A);
    A = B;
    rank += npiv;
    n = A->n;
    m = A->m;

    npiv = spasm_find_pivots(A, p, qinv);
    density = spasm_schur_probe_density(A, p, qinv, npiv, 100);
  }

  /* sparse schur complement : GPLU */
  if (density < 0.1) {
    spasm_lu *LU = spasm_LU(A, p, SPASM_DISCARD_L);
    rank +=  LU->U->n;
    spasm_free_LU(LU);
  } else {
    /* dense schur complement */
    int r = spasm_schur_rank(A, p, npiv);
    fprintf(stderr, "rank = %d + %d\n", npiv, r);
    rank += npiv + r;
  }

  end_time = spasm_wtime();
  fprintf(stderr, "done in %.3f s rank = %d\n", end_time - start_time, rank);
  spasm_csr_free(A);
  free(p);
  free(qinv);
  return 0;
}