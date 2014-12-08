#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm_triplet *T;
  spasm *C;
  int i, n, test, u, v;
  int *p, *q;

  assert(argc == 2);
  test = atoi(argv[1]);

  T = spasm_load_sms(stdin, 42013);
  C = spasm_compress(T);
  spasm_triplet_free(T);
  n = C->n;

  p = spasm_row_sort(C);
  q = spasm_malloc(n * sizeof(int));

  for(i = 0; i < n; i++) {
    q[i] = 0;
  }
  for(i = 0; i < n; i++) {
    q[ p[i] ]++;
    printf("%d\t%d\t%d\n", i, p[i], spasm_row_weight(C, p[i]));
  }

  // test A
  for(i = 0; i < n; i++) {
    if (q[i] != 1) {
      printf("not ok %d - sort not bijective\n", test);
      exit(1);
    }
  }

  // test B
  v = spasm_row_weight(C, p[0]);
  for(i = 1; i < n; i++) {
    u = spasm_row_weight(C, p[i]);
    if (u < v) {
      printf("not ok %d - sort not increasing\n", test);
      exit(1);
    }
    v = u;
  }

  printf("ok %d - sort\n", test);

  spasm_csr_free(C);
  free(p);
  free(q);
  return 0;
}
