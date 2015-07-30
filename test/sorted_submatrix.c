#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *S;
  int n, m, k;
  int *py, *Ap;

  T = spasm_load_sms(stdin, 46337);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;
  Ap = A->p;
  py = spasm_malloc(n * sizeof(int));

  spasm_row_entries_sort(A, 1);
  for(k = 0; k < n; k++) {
    py[k] = Ap[k];
  }

  assert(n > 10);
  assert(m > 10);
  S = sorted_spasm_submatrix(A, 0, 5, 0, 5, py, SPASM_WITH_NUMERICAL_VALUES);
  spasm_save_csr(stdout, S);
  spasm_csr_free(S);

  S = sorted_spasm_submatrix(A, 0, 5, 5, 10, py, SPASM_WITH_NUMERICAL_VALUES);
  spasm_save_csr(stdout, S);
  spasm_csr_free(S);
  free(py);
  spasm_csr_free(A);
}
