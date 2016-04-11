#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main(){
  
  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, i;
  spasm_triplet * T = spasm_load_sms(stdin, prime);
  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

  /* find free pivots. */
  int n_cheap;
  int *p = spasm_cheap_pivots(A, &n_cheap);

  /* filter rows that are too dense */
  int free_nnz = 0;
  int min = A->m;
  int max = 0;
  for(i=0; i<n_cheap; i++) {
      free_nnz += spasm_row_weight(A, p[i]);
      min = spasm_min(min, spasm_row_weight(A, p[i]));
      max = spasm_max(max, spasm_row_weight(A, p[i]));
  }
  float avg = 1.0 * free_nnz / n_cheap;

  fprintf(stderr, "[schur] free pivots NNZ (min/avg/max): %d / %.1f / %d\n", min, avg, max);
  
  for(i=0; i<n_cheap; i++) {
    int inew = p[i];

    if (spasm_row_weight(A, inew) > 3*avg) {
      /* this row is too dense */
      spasm_swap(p, inew, --n_cheap);
      i--;
    } 
  }
  
  fprintf(stderr, "[schur] %d free pivots after filtering\n", n_cheap);


  /* compute schur complement */
  spasm *S = spasm_schur(A, p, n_cheap);

  fprintf(stderr, "Schur complement: (%d x %d), nnz : %d, dens : %.5f\n", S->n, S->m, spasm_nnz(S), 1. * spasm_nnz(S)/(1.*S->n * S->m));

  spasm_save_csr(stdout, S);

  spasm_csr_free(S);
  spasm_csr_free(A);
}
