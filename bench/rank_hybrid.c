#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/** Computes the rank of the input matrix using the hybrid strategy */

/* NOT DRY (the same function is in rank_hybrid) */
spasm * filtered_schur(spasm *A, int *npiv){
  int n_cheap, n_filtered, free_nnz, min, max, h, i;
  float avg;

 /* find free pivots. */
  int *p = spasm_cheap_pivots(A, &n_cheap);
  int *filtered = malloc(A->n * sizeof(int));

  /* collect stats */
  free_nnz = 0;
  min = A->m;
  max = 0;
  for(i=0; i<n_cheap; i++) {
    h =  spasm_row_weight(A, p[i]);
    free_nnz += h;
    min = spasm_min(min, h);
    max = spasm_max(max, h);
  }
  avg = 1.0 * free_nnz / n_cheap;
  fprintf(stderr, "[schur] free pivots NNZ (min/avg/max): %d / %.1f / %d\n", min, avg, max);

  /* filter rows that are too dense */
  n_filtered = 0;
  h = n_cheap-1;
  for(i=0; i<n_cheap; i++) {
    if (spasm_row_weight(A, p[i]) <= 3*avg) {
      filtered[n_filtered++] = p[i];
    } else {
      filtered[h--] = p[i];
    }
  }
  for(i=n_cheap; i<A->n; i++) {
    filtered[i] = p[i];
  }
  
  fprintf(stderr, "[schur] %d free pivots after filtering\n", n_filtered);
  free(p);

  /* schur complement */
  spasm *S = spasm_schur(A, filtered, n_filtered);

  fprintf(stderr, "Schur complement: (%d x %d), nnz : %d, dens : %.5f\n", S->n, S->m, spasm_nnz(S), 1. * spasm_nnz(S)/(1.*S->n * S->m));

  free(filtered);
  *npiv = n_filtered;

  return S;
}


int main(int argc, char **argv) {
  
  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, n_times, i, n_cheap, rank, npiv;
  double start_time, end_time;
  spasm_triplet * T;
  spasm *A, *B;

  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  /* 3 iterations of Schur complement, by default */
  n_times = 3;
  if(argc > 1){
    n_times = atoi(argv[1]);
  }

  start_time = spasm_wtime();

  rank = 0;
  for(i = 0; i < n_times; i++){
    fprintf(stderr, "%d : A : (%d x %d) nnz %d\n", i, A->n, A->m, A->nzmax);
    B = filtered_schur(A, &npiv);
    spasm_csr_free(A);
    A = B;
    rank += npiv;
  }

  /* finish the job with GPLU */
  int *p = spasm_cheap_pivots(A, &n_cheap);
  spasm_lu *LU = spasm_LU(A, p, 0);
  free(p);

  end_time = spasm_wtime();
  fprintf(stderr,"done in %.3f s rank = %d\n", end_time - start_time, rank + LU->U->n);

  spasm_free_LU(LU);
  spasm_csr_free(A);
}
