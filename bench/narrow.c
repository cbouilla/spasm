#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main(){
  
  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, i, n_cheap, n_filtered, free_nnz, min, max, h;
  float avg;
  spasm_triplet * T = spasm_load_sms(stdin, prime);
  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

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

  /* narrow schur trick */

 int r = spasm_narrow_schur_trick(A, filtered, n_filtered);
  printf("rank : %d\n", r);

  free(filtered);
  spasm_csr_free(A);


}
