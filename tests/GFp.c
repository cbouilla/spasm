#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main() {
  spasm_GFp prime = 42013;
  
  for (spasm_GFp i = 1; i < prime; i++) {
    spasm_GFp j = spasm_GFp_inverse(i, prime);
    spasm_GFp k = (i * j) % prime;
    fprintf(stderr, "1/%d == %d mod %d\n", i, j, prime);
    assert(k == 1);
  }
  exit(EXIT_SUCCESS);
}
