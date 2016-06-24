/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_column.c  */
#include <assert.h>
#include "spasm.h"

/* Given a matrix A, permute the columns of A such that the leftmost columns are
 * sparsest than the right ones */

spasm * spasm_sort_columns(spasm *A, int **col_perm_pt){
  spasm *B;
  int i, px, m, n, cnz_max, *Ap, *Aj, *Hw, *w, *col_perm;

  // check inputs
  assert(A != NULL);

  m = A->m;
  n = A->n;
  Ap = A->p;
  Aj = A->j;
  
  /*get and initialize workspace*/
  col_perm = spasm_malloc(m * sizeof(int)); // permutation vector
  Hw = spasm_malloc(m * sizeof(int)); // Hamming weight of each columns
  *col_perm_pt = col_perm;
  
  for(i = 0; i < m; i++){
    col_perm[i] = 0;
    Hw[i] = 0;
  }
  cnz_max = 0; //number of non zero entries in the densest column

  /* Browse matrix A to cumpute the Hamming weight of each columns */
  for(i = 0; i < n; i++){
    for(px = Ap[i]; px < Ap[i+1]; px++){
      Hw[Aj[px]]++;
      if(Hw[Aj[px]] > cnz_max) cnz_max = Hw[Aj[px]];
    }
  }

  assert(cnz_max <= n);

  //get workspace
  w = spasm_malloc((cnz_max + 1) * sizeof(int)); // pointer on column index.
  
  for(i = 0; i < cnz_max + 1; i++){
    w[i] = 0;
  }

  for(i = 0; i < m; i++){
    if(Hw[i] != cnz_max){
      w[Hw[i]+1]++; // count the number of columns that have Hw[i] entries.
    }
  }

  for(i = 0; i < cnz_max; i++){
    w[i+1] += w[i]; // w[i] <--- first column that has exactly Hw[i] entries.
  }

  /* write permutation vector */
  for(i = 0; i < m; i++){
    col_perm[i] = w[Hw[i]];
    w[Hw[i]]++; // increment pointer.
  }

  // free workspace
  free(w);
  free(Hw);

  // permute matrix.
  B = spasm_permute(A, NULL, col_perm, (A->x != NULL));

  return B;
}

