#include <assert.h>
#include <stdlib.h>

#include "spasm.h"

/*
 * Solve x.A = b
 * 
 * b has size m (#columns of A), solution has size n.
 * 
 * returns true if a solution exists
 */
bool spasm_solve(const spasm_lu *fact, const spasm_ZZp *b, spasm_ZZp *x)
{
	const spasm *L = fact->L;
	const spasm *U = fact->U;
	assert(L != NULL);
	// int n = L->n;
	int m = U->m;
	int r = U->n;   /* rank */

	/* get workspace */
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	spasm_ZZp *z = spasm_malloc(r * sizeof(*z));
	
	/* inverse permutation for U */
	int *Uq = spasm_malloc(r * sizeof(*Uq));
	const int *Uqinv = fact->Uqinv;
	for (int j = 0; j < m; j++) {
		int i = Uqinv[j];
		if (i != -1)
			Uq[i] = j;
	}

	/* z.U = b  (if possible) */
	for (int i = 0; i < m; i++)
		y[i] = b[i];
	bool ok = spasm_dense_forward_solve(U, y, z, Uq);

	/* y.LU = b */
	spasm_dense_back_solve(L, z, x, fact->Lqinv);
	
	free(y);
	free(z);
	free(Uq);
	return ok;
}


spasm * spasm_solve_gesv(const spasm_lu *fact, const spasm *b)
{
	(void) fact;
	(void) b;
	return NULL;
}
