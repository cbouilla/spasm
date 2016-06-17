#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/** Finds "cheap pivots" (i.e. Faugère-Lachartre pivots) and computes the Schur complement w.r.t. these pivots */

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

int main(){
  int n_piv, prime = 42013;
  spasm_triplet *T; 
  spasm *A, *S;

  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  S = filtered_schur(A, &n_piv);
  spasm_save_csr(stdout, S);
  spasm_csr_free(S);
  spasm_csr_free(A);
  return 0;
}
