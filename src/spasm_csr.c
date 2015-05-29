#include <assert.h>
#include <stdbool.h>
#include "spasm.h"

void spasm_row_entries_sort(spasm *M, int with_value) {
  int i, k, l, n, *p, *j, tmpj;
  spasm_GFp *x, tmpx;

  n = M->n; // <--- number of rows
  p = M->p; // <--- row pointers
  j = M->j; // <--- column index
  x = M->x; // <--- matrix entries

  for(i = 0; i < n; i++) {
    for(k = p[i] + 1; k < p[i+1] ; k++) {
      tmpj = j[k];
      if (with_value) tmpx = x[k];
      l = k-1;
      while (l >= 0 && j[l] > tmpj) {
	j[l+1] = j[l];
	if (with_value) x[l+1] = x[l];
	l--;
      }
      j[l+1] = tmpj;
      if (with_value) x[l+1] = tmpx;
    }
  }

}
