#include <stdio.h>
#include <stdlib.h>

#include "spasm.h"

int main() {
  int prime, fail;
  spasm_GFp i, j;

  fail = 0;
  prime = 257;
  for(i = 1; i < prime; i++) {
    j = spasm_GFp_inverse(i, prime);
    fail |= (((j * i) % prime) != 1);
  }
  if (fail) {
    printf("not ok - GFp inverse\n");
    exit(1);
  } else {
    printf("ok - GFp inverse\n");
  }
  exit(EXIT_SUCCESS);
}
