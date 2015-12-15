#include<assert.h>
#include"spasm.h"

/*
 * Given a diagonal index d>0, index k and i, 
 * permutations p and R, pivots numbers indicators 
 * ri and rj, return the corresponding index k for
 * the next stage (solving y_k L_{d-1,k} = x_k)
 * and the position of i in vector x_k_new.
 */
int spasm_lazy_solve_update(int d, int k, int n_blocks, int *i_ptr, const int **ri, const int **rj, const int **p){
  int i, d_new, dr, bound; 

  //check inputs.
  assert(d > 0);
  assert(k >= 0);
  assert(k+1 < n_blocks-d);

  i = *i_ptr;
  bound = rj[k+d][d-1];

  if(d > 1){
    if(i < bound) {
      k++;
    }
    else {
      d_new = d-1;
      i = i - bound;
      dr = ri[k][d_new] - ri[k][d_new];
      i = i + rj[k+d][d_new-1];
      i = i + dr;
      *i_ptr = i;
    }
  }
  else {
    if (i < bound){
      k++;
      i = p[k][i];
      *i_ptr = i;
    }
    else {
      i = i - bound;
      i = i + ri[k][0];
      i = p[k][i];
      *i_ptr = i;
    }
  }
  return k;
}
