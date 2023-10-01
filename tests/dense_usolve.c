#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv)
{
  struct spasm_triplet *T = spasm_triplet_load(stdin, 32003, NULL);
  struct spasm_csr *G = spasm_compress(T);
  spasm_triplet_free(T);

  int n = G->n;
  int m = G->m;
  assert(n <= m);

  spasm_ZZp *x = malloc(n * sizeof(spasm_ZZp));
  spasm_ZZp *b = malloc(m * sizeof(spasm_ZZp));
  spasm_ZZp *y = malloc(m * sizeof(spasm_ZZp));

  /* test A ------------------------- with a sensible RHS ----------- */
  printf("# testing correct solution\n");

  for (int i = 0; i < n; i++)
    x[i] = spasm_ZZp_init(G->field, rand());
  
  for (int j = 0; j < m; j++)
    b[j] = 0;
  
  spasm_xApy(x, G, b);
  for (int j = 0; j < m; j++)
    y[j] = b[j];

  bool result = spasm_dense_forward_solve(G, y, x, SPASM_IDENTITY_PERMUTATION);
  if (!result) {
    printf("not ok - dense forward-substitution triangular solver [solution not found]\n");
    exit(1);
  }

  for (int j = 0; j < m; j++)
    y[j] = 0;
 
  spasm_xApy(x, G, y);
  for (int i = 0; i < m; i++)
    if (y[i] != b[i]) {
      printf("not ok - dense forward-substitution triangular solver [incorrect solution found]\n");
      exit(1);
    }
  

  /* test B ------------------------- with a bogus RHS ----------- */
  if (n < m) {
    printf("# testing bogus solution\n");
    for (int j = 0; j < m; j++)
      b[j] = spasm_ZZp_init(G->field, rand());

    result = spasm_dense_forward_solve(G, b, x, SPASM_IDENTITY_PERMUTATION);
    if (result) {
      printf("not ok - dense forward-substitution triangular solver [bogus solution found]\n");
      exit(1);
    }
  }

  printf("ok -  dense forward-substitution triangular solver\n");

  spasm_csr_free(G);
  free(x);
  free(y);
  free(b);
  exit(EXIT_SUCCESS);
}
