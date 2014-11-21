#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
  spasm *T, *G;
  int i, n, root, *xi, *pstack;

  assert(argc > 1);
  root = atoi(argv[1]);

  T = spasm_load_ctf(stdin, 257);
  G = spasm_compress(T);
  spasm_spfree(T);

  n = G->n;

  pstack = malloc(n * sizeof(int));
  xi = malloc(n * sizeof(int));

  i = spasm_dfs(root, G, n, xi, pstack, NULL);
  for( ; i < n; i++) {
    printf("%d\n", xi[i]);
  }

  return 0;
}
