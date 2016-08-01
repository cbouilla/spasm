/* indent -nfbs -i2 -nip -npsl -di0 -nut structural_rank.c  */
#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/* Provide a lower-bound on the rank of the input matrix without performing any arithmetic operation */

int main(int argc, char **argv) {
  /* charge la matrice depuis l'entrée standard */
  int npiv, *p, *qinv, n, m, *Aj, nnz;

  assert(argc > 1);
  spasm *A = spasm_load_CADO(argv[1]);
  Aj = A->j;
  m = A->m;
  n = A->n;
  nnz = A->nzmax;

  #pragma omp parallel for schedule(static)
  for(int i=0; i< nnz; i++)
    Aj[i] = (m-1) - Aj[i];

  p = spasm_malloc(n * sizeof(int));
  qinv = spasm_malloc(m * sizeof(int));
  npiv = spasm_find_pivots(A, p, qinv);

  free(p);
  free(qinv);
  spasm_csr_free(A);
}
