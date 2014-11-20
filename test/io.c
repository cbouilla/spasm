#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm *T, *C;
  int test;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_ctf(stdin, 257);
  switch(test) {
  case 1:
    spasm_save_ctf(stdout, T);
    spasm_spfree(T);
    break;

  case 2:
    C = spasm_compress(T);
    spasm_spfree(T);
    spasm_save_ctf(stdout, C);
    spasm_spfree(C);
  }
}
