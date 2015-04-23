#include <assert.h>
#include "spasm.h"


/* x = x + beta * A[j], where x is a dense vector and A[j] is sparse
 *
 * low-level operation for maximum flexibility;
 *
 * This is where all the heavy lifting should take place.
 */
void spasm_scatter(const int *Aj, const spasm_GFp *Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime) {
  int j, p;

    for (p = from; p < to; p++) {
      j = Aj[p];
      // axpy-inplace
      x[j] = (x[j] + ((beta * Ax[p]))) % prime /* ultra-naive */;
    }
    
}
