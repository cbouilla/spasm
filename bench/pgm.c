#include <stdio.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *A;
  int x, y;

  if (argc < 3) {
    printf("USAGE: %s xdim ydim < matrix.sms > file.pgm\n", argv[0]);
    exit(1);
  }

  x = atoi(argv[1]);
  y = atoi(argv[2]);

  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  spasm_save_pgm(stdout, x, y, A);

  spasm_csr_free(A);
  return 0;
}
