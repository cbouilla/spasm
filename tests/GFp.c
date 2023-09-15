#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main() {
  i64 prime = 42013;
  spasm_field F;
  spasm_field_init(prime, &F);
  
  for (i64 i = 1; i < prime; i++) {
    spasm_ZZp x = spasm_ZZp_init(&F, i);
    spasm_ZZp j = spasm_ZZp_inverse(&F, x);
    spasm_ZZp k = spasm_ZZp_mul(&F, i, j);
    fprintf(stderr, "%" PRId64 " * %d mod == %d mod %" PRId64 "\n", i, j, k, prime);
    assert(k == 1);
    assert(j <= prime / 2);
    assert(j >= -prime / 2);
  }
  exit(EXIT_SUCCESS);
}
