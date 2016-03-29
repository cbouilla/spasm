#include<assert.h>
#include"spasm.h"

/*
 * Given a matrix A super_spasm, a list of d matrices L super_sparse,
 * a row i of A (in big matrix), solve the successive "lazy" system x * L = y.
 * (If d = 0, L is NULL). Compute the product u = x * A.
 *
 * The return value unz is the number of non-zero entries in U.
 */
int super_spasm_lazy(super_spasm *A, super_list *L, int i, spasm_GFp *u, int *ui){
  int n_big, top, end, j, k, tmpi, tmpy, unz;
  int *xi, *yi, *Aperm;
  spasm_GFp *x, *y;
  super_list *Ltmp;

  /*check inputs */
  assert(A != NULL);
  n_big = A->n;
  Ltmp = L;
  if(Ltmp != NULL){
    assert(Ltmp->S->n = n_big); // check size compatibility.
  }
 
  Aperm = A->p;

  /* get workspace */
  x = (Ltmp != NULL) ? spasm_malloc(n_big * sizeof(spasm_GFp)) : NULL;
  xi = (Ltmp != NULL) ? spasm_malloc(3 * n_big * sizeof(int)) : NULL;
  y = spasm_malloc(n_big * sizeof(spasm_GFp));
  yi = spasm_malloc(n_big * sizeof(int));

  /* initialize workspace */

  //clear workspace.
  for(j = 0; j < n_big; j++){
    if(Ltmp != NULL){
      x[j] = 0;
    }
    y[j] = 0;
    yi[j] = 0;
  }
  if(Ltmp != NULL){
    for(j = 0; j < 3 * n_big; j++){
      xi[j] = 0;
    }
  }

  end = 1; // yi pattern ends at 1.
  yi[0] = i;
  y[yi[0]] = 1; // y = e_i.

  /* main loop : solve successive system */
  while(Ltmp != NULL){
    //solve system x * L = y 
    top = super_spasm_sparse_solve(Ltmp->S, y, yi, end, x, xi);
    
    //clear workspace :
    for(j = 0; j < n_big; j++){
      y[j] = 0;
      yi[j] = 0;
    }
    //scatter x in y.
    for(j = top; j < n_big; j++){
      yi[j - top] = xi[j];
      y[yi[j - top]] = x[xi[j]];
      //clear x and xi :
      x[xi[j]] = 0;
      xi[j] = 0;
    }
    end = n_big - top;
    //update L
    Ltmp = Ltmp->next;
  }

  /* compute the product u = y * A */
  // first sort entries in y : insertion sort.
  for(j = 1; j < end; j++){
    for(k = j; k > 0 && yi[k] < yi[k - 1]; k--){
      tmpi = yi[k];
      tmpy = y[yi[k]];
      yi[k] = yi[k - 1];
      y[yi[k]] = y[yi[k - 1]];
      yi[k - 1] = tmpi;
      y[yi[k - 1]] = tmpy;
    }
  }

  unz = super_sparse_gax(A, y, yi, end, u, ui);

  /* free workspace */
  free(y);
  free(yi);
  if(xi != NULL) free(xi);
  if(x != NULL) free(x);

  return unz;
}
