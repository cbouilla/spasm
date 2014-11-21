#include <assert.h>
#include "spasm.h"


/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse
 *
 * low-level operation for maximum flexibility;
 *
 * This is where all the heavy lifting should take place.
 */
void spasm_scatter(const int *Ai, const spasm_GFp *Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime) {
  int i, p;

    for (p = from; p < to; p++) {
      i = Ai[p];
      // axpy-inplace
      x[i] = (x[i] + ((beta * Ax[p]) % prime)) % prime /* ultra-naive */;
    }
}
