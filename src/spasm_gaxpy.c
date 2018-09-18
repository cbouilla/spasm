#include <assert.h>
#include "spasm.h"

/*
 * (dense vector) * (sparse) Matrix y <--- y + x*A
 */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp * y) {
	const int n = A->n;
	const int *Ap = A->p;
	const int *Aj = A->j;
	const spasm_GFp *Ax = A->x;
	const int prime = A->prime;
	assert(Ax != NULL);

	for (int i = 0; i < n; i++)
		spasm_scatter(Aj, Ax, Ap[i], Ap[i + 1], x[i], y, prime);
}

