#include <assert.h>
#include <stdio.h>
#include "spasm.h"


int main(){  
  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, n_cheap, *p;
  double start_time, end_time;

  spasm_triplet * T = spasm_load_sms(stdin, prime);
  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

  start_time = spasm_wtime();
  p = spasm_cheap_pivots(A, &n_cheap);
  end_time = spasm_wtime();
  printf("%d; %.2f\n", n_cheap, end_time - start_time);
  
  free(p);
  spasm_csr_free(A);
}
