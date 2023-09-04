#include "spasm.h"
/*
 * x = x + beta * A[j], where x is a dense vector and A[j] is sparse
 * 
 * low-level operation for maximum flexibility;
 * 
 * This is where all the heavy lifting should take place.
 */
void spasm_scatter(const int *Aj, const spasm_GFp * Ax, i64 from, i64 to, spasm_GFp beta, spasm_GFp * x, int prime) {
	for (i64 px = from; px < to; px++) {
		int j = Aj[px];
		x[j] = (x[j] + ((beta * Ax[px]))) % prime; /* ultra-naive */ 
	}
}
