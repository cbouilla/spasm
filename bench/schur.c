#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main(){
  
  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, i, n_cheap, n_filtered, free_nnz, min, max;
  float avg;
  spasm_triplet * T = spasm_load_sms(stdin, prime);
  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

  /* find free pivots. */
  int *p = spasm_cheap_pivots(A, &n_cheap);
  int *filtered = malloc(n_cheap * sizeof(int));

  /* collect stats */
  free_nnz = 0;
  min = A->m;
  max = 0;
  for(i=0; i<n_cheap; i++) {
      free_nnz += spasm_row_weight(A, p[i]);
      min = spasm_min(min, spasm_row_weight(A, p[i]));
      max = spasm_max(max, spasm_row_weight(A, p[i]));
  }
  avg = 1.0 * free_nnz / n_cheap;
  fprintf(stderr, "[schur] free pivots NNZ (min/avg/max): %d / %.1f / %d\n", min, avg, max);

  /* filter rows that are too dense */
  n_filtered = 0;
  for(i=0; i<n_cheap; i++) {
    if (spasm_row_weight(A, p[i]) <= 3*avg) {
      filtered[n_filtered++] = p[i];
    }
  }
  
  fprintf(stderr, "[schur] %d free pivots after filtering\n", n_filtered);
  free(p);

  /* compute schur complement */
  spasm *S = spasm_schur(A, filtered, n_filtered);

  fprintf(stderr, "Schur complement: (%d x %d), nnz : %d, dens : %.5f\n", S->n, S->m, spasm_nnz(S), 1. * spasm_nnz(S)/(1.*S->n * S->m));

  free(filtered);
  spasm_csr_free(A);

  spasm_save_csr(stdout, S);
  spasm_csr_free(S);
}
