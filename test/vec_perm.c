#include <stdio.h>
#include "spasm.h"

static inline void swap(int *a, int i, int j) {
  int x = a[i];
    a[i] = a[j];
    a[j] = x;
}

int main(int argc, char **argv) {
  int n = 50;
  int i;

  printf("1..2\n");

  int *p = malloc(n * sizeof(int));
  int *pinv;
  spasm_GFp *x = malloc(n * sizeof(spasm_GFp));
  spasm_GFp *y = malloc(n * sizeof(spasm_GFp));

  for(i = 0; i < n; i++) {
    x[i] = i*i + 3*i - 7;
    p[i] = i;
  }

  // generate random permutation
  for(i = n-1; i > 0; i--) {
    swap(p, i, rand() % i);
  }

  // test 1 : apply permutation
  cs_pvec(p, x, y, n);
  int fail = 0;
  for(i = 0; i < n; i++) {
    fail |= (y[i] == i*i + 3*i - 7);
  }
  if (fail) {
    printf("not ok 1 - vector permutation\n");
  } else {
    printf("ok 1 - vector permutation\n");
  }

  // test 2 : apply inverse permutation
  pinv = spasm_pinv(p, n);
  cs_pvec(pinv, y, x, n);
  for(i = 0; i < n; i++) {
    fail |= (x[i] != i*i + 3*i - 7);
  }
  if (fail) {
    printf("not ok 2 - inverse vector permutation\n");
  } else {
    printf("ok 2 - inverse vector permutation\n");
  }

  // todo : test matrix permutation

}
