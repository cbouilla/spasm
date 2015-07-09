#include <stdio.h>
#include <assert.h>
#include "spasm.h"


int main() {
  spasm_triplet *T;
  spasm *A, *B, *C;
  int m, k;
 
 
  // Load matrix A
  T = spasm_load_sms(stdin, 46337);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  m = A->m;
 
  // load B
  T = spasm_triplet_alloc(5, m, 10, 46337, 1);
  spasm_add_entry(T, 0, 0, 1);
  spasm_add_entry(T, 0, m/2, 2);
  for (k = 0; k < 5; k++) {
    spasm_add_entry(T, k, m-1-k, k+1);
  }
  B = spasm_compress(T);
  spasm_triplet_free(T);

  C = spasm_row_concatenation(A, B, 1);
  spasm_save_csr(stdout, C);

  spasm_csr_free(C);
  spasm_csr_free(A);
  spasm_csr_free(B);

  return 0;
}
