/* indent -nfbs -i2 -nip -npsl -di0 -nut structural_rank.c  */
#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/* Provide a lower-bound on the rank of the input matrix without performing any arithmetic operation */

int main(int argc, char **argv) {
  /* charge la matrice depuis l'entrée standard */
  int n_cheap, *p;
  double start_time, end_time;

  assert(argc > 1);
  spasm *A = spasm_load_CADO(argv[1]);
  
  for(int i=0; i<A->nzmax; i++) {
    A->j[i] = A->m - A->j[i];
  }

  start_time = spasm_wtime();
  p = spasm_cheap_pivots(A, &n_cheap);
  end_time = spasm_wtime();
  // printf("%d; %.2f\n", n_cheap, end_time - start_time);

  spasm_save_csr(stdout, A);

  free(p);
  spasm_csr_free(A);
}
