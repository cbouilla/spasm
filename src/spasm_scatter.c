#include "spasm.h"
/*
 * x = x + beta * A[i], where x is a dense vector
 * 
 * This is where all the heavy lifting should take place.
 */
void spasm_scatter(const spasm *A, int i, spasm_GFp beta, spasm_GFp * x)
{
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_GFp *Ax = A->x;
	const spasm_GFp prime = A->prime;
	for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
		int j = Aj[px];
		x[j] = (x[j] + beta * Ax[px]) % prime; /* ultra-naive */ 
	}
}
