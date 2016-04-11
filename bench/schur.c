#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main(){
  
  // charge la matrice depuis l'entrée standard
  int prime = 42013;
  spasm_triplet * T = spasm_load_sms(stdin, prime);
  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

  // permutation p : cheap pivots.
  int n_cheap;
  int *p = spasm_cheap_pivots(A, &n_cheap);

  // calculer le complément de schur après n_cheap.
  spasm *S = spasm_schur(A, p, n_cheap);
  fprintf(stderr, "Schur : (%d x %d), nnz : %d, dens : %.5f\n", S->n, S->m, spasm_nnz(S), 1. * spasm_nnz(S)/(S->n * S->m));

  spasm_save_csr(stdout, S);

  spasm_csr_free(S);
  spasm_csr_free(A);
}
