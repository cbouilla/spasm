/* indent -nfbs -i2 -nip -npsl -di0 -nut spasm_scatter.c  */
#include <assert.h>
#include "spasm.h"

/* given a permutation vector col_perm, permute the columns of matrix A */
spasm * spasm_columns_permute(spasm *A, int *col_perm){
  spasm *B;
  int n, m, nnz, prime, *Ap, *Aj, with_values, i, px, *Bp, *Bj;
  spasm_GFp *Ax, *Bx;

  /*check inputs*/
  assert(A != NULL);
  if(col_perm == NULL) return A; // no permutation.

  n = A->n;
  m = A->m;
  nnz = A->nzmax;
  prime = A->prime;
  Ap = A->p;
  Aj = A->j;
  with_values = (A->x != NULL) ? 1 : 0;

  /* allocate result*/
  B = spasm_csr_alloc(n, m, nnz, prime, with_values);
  Bp = B->p;
  Bj = B->j;
  Bx = (with_values) ? B->x : NULL;  

  /* write output matrix */
  for(i = 0; i < n; i++){
    Bp[i] = Ap[i];
    for(px = Ap[i]; px < Ap[i+1]; px++){
      Bj[px] = col_perm[Aj[px]]; // new column index.
      if(with_values){
	Bx[px] = A->x[px];
      }
    }
  }

  /* finalize and return B */
  Bp[n] = Ap[n];
  return B;
}

/* Given a matrix A, permute the columns of A such that the leftmost columns are
 * sparsest than the right ones */

spasm * spasm_sort_columns(spasm *A){
  spasm *B;
  int i, px, m, n, cnz_max, *Ap, *Aj, j, *Hw, *w, *col_perm;

  // check inputs
  assert(A != NULL);

  m = A->m;
  n = A->n;
  Ap = A->p;
  Aj = A->j;
  
  /*get and initialize workspace*/
  col_perm = spasm_malloc(m * sizeof(int)); // permutation vector
  Hw = spasm_malloc(m * sizeof(int)); // Hamming weight of each columns

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
  B = spasm_columns_permute(A, col_perm);
  free(col_perm);

  return B;
}
