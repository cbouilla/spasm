#include <assert.h>
#include <stdio.h>
#include "spasm.h"


int main(){  
  /* charge la matrice depuis l'entrée standard */
  int prime = 42013, n_cheap, n, m, i, j, k, top;
  int *p, *q, *qinv, *Ap, *Aj, *marks, *xj, *pstack;
  spasm *A, *B;
  spasm_triplet *T;

  T = spasm_load_sms(stdin, prime);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

#if 0
  /* DM preprocessing */
  spasm_dm * x = spasm_dulmage_mendelsohn(A);
  qinv = spasm_pinv(x->DM->q, A->m);
  B = spasm_permute(A, x->DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  free(qinv);
  spasm_csr_free(A);
  A = B;
#endif

  Ap = A->p;
  Aj = A->j;

  /* start_time = spasm_wtime(); */
  p = spasm_cheap_pivots(A, &n_cheap);
  /* end_time = spasm_wtime(); */

  /* build qinv to reflect the changes in p */
  qinv = spasm_malloc(m * sizeof(int));
  for(j = 0; j < m; j++) {
    qinv[j] = -1;
  }
  
  /* pivotal column first, in row-order */
  k = 0;
  for(i = 0; i < n_cheap; i++) {
    j = Aj[ Ap[ p[i] ] ];   /* the pivot is the first entry of each row */
    assert(qinv[j] == -1);
    qinv[j] = k++;
  }

  /* put remaining non-pivotal columns afterwards, in any order */
  for(j = 0; j < m; j++) {
    if (qinv[j] == -1) {
      qinv[j] = k++;
    }
  }
  // printf("m = %d / k = %d\n", m, k);
  assert(k == m);

  B = spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  spasm_save_csr(stdout, B);

  free(p);
  // free(q);
  free(qinv);
  spasm_csr_free(A);
  spasm_csr_free(B);

  return 0;
}
