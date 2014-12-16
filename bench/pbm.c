#include <stdio.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A;


  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  spasm_save_pbm(stdout, A);

  spasm_csr_free(A);
  return 0;
}
