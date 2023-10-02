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
bool spasm_solve(const struct spasm_lu *fact, const spasm_ZZp *b, spasm_ZZp *x)
{
	const struct spasm_csr *L = fact->L;
	const struct spasm_csr *U = fact->U;
	assert(L != NULL);
	// int n = L->n;
	int m = U->m;
	int r = U->n;   /* rank */

	/* get workspace */
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	spasm_ZZp *z = spasm_malloc(r * sizeof(*z));
	
	/* inverse permutation for U */
	int *Uq = spasm_malloc(r * sizeof(*Uq));
	const int *qinv = fact->qinv;
	for (int j = 0; j < m; j++) {
		int i = qinv[j];
		if (i != -1)
			Uq[i] = j;
	}

	/* z.U = b  (if possible) */
	for (int i = 0; i < m; i++)
		y[i] = b[i];
	bool ok = spasm_dense_forward_solve(U, y, z, Uq);

	/* y.LU = b */
	spasm_dense_back_solve(L, z, x, fact->p);
	
	free(y);
	free(z);
	free(Uq);
	return ok;
}

/* Solve XA == B (returns garbage if a solution does not exist).
 * If ok != NULL, then sets ok[i] == 1 iff xA == B[i] has a solution
 */
struct spasm_csr * spasm_gesv(const struct spasm_lu *fact, const struct spasm_csr *B, bool *ok)
{
	i64 prime = B->field->p;
	assert(prime == fact->L->field->p);
	assert(fact->L != NULL);
	int n = B->n;
	int m = B->m;
	int Xm = fact->L->n;
	struct spasm_triplet *X = spasm_triplet_alloc(n, Xm, (i64) Xm * n, prime, true);
	int *Xi = X->i;
	int *Xj = X->j;
	spasm_ZZp *Xx = X->x;

	#pragma omp parallel
	{
		spasm_ZZp *b = spasm_malloc(m * sizeof(*b));
		spasm_ZZp *x = spasm_malloc(Xm * sizeof(*x));
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) 
				b[j] = 0;
			spasm_scatter(B, i, 1, b);
			bool res = spasm_solve(fact, b, x);
			if (ok)
				ok[i] = res;
			for (int j = 0; j < Xm; j++) 
				if (x[j] != 0) {
					i64 xnz;
					#pragma omp atomic capture
					{ xnz = X->nz; X->nz += 1; }
					Xi[xnz] = i;
					Xj[xnz] = j;
					Xx[xnz] = x[j];
				}
		}
		free(b);
		free(x);
	}
	struct spasm_csr *XX = spasm_compress(X);
	spasm_triplet_free(X);
	return XX;
}
