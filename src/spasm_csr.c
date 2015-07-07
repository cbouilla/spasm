#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

void spasm_row_entries_sort(spasm *M, int with_value) {
  int i, k, l, n, *Mp, *Mj, tmpj;
  spasm_GFp *Mx, tmpx;

  n = M->n; // <--- number of rows
  Mp = M->p; // <--- row pointers
  Mj = M->j; // <--- column index
  Mx = M->x; // <--- matrix entries

  for(i = 0; i < n; i++) {
    for(k = Mp[i] + 1; k < Mp[i+1] ; k++) {
      tmpj = Mj[k];
      if (with_value) tmpx = Mx[k];
      for(l = k; l > Mp[i] && Mj[l-1] > tmpj; l--) {
	Mj[l] = Mj[l-1];
	if (with_value) Mx[l] = Mx[l-1];
      }
      Mj[l] = tmpj;
      if (with_value) Mx[l] = tmpx;
    }
  }

  
}
