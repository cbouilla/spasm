#include <assert.h>
#include "spasm.h"

/*
 * Returns a matrix K whose rows span a basis of the RIGHT kernel of A, i.e. a basis of the vector
 * space {x | A.x == 0 }.
 *
 * So, in principle K . At == 0.
 */
spasm *spasm_kernel(spasm * A) 
{
	assert(A->x);

	spasm_lu *LU = spasm_echelonize(A, -1); // careful, this destroys *A
	spasm *U = LU->U;
	spasm *L = spasm_transpose(U, SPASM_WITH_NUMERICAL_VALUES);

	int n = L->n;
	int m = L->m;
	int *qinv = LU->qinv;
	int prime = L->prime;

	/* in L, pivots must be pushed to the first entry of pivotal rows for the sparse triangular solver */
	int *Lp = L->p;
	int *Lj = L->j;
	spasm_GFp *Lx = L->x;
	int *q = spasm_malloc(m * sizeof(*q));

	for (int i = 0; i < m; i++) {
		if (qinv[i] < 0)   /* non-pivotal row */
			continue;
		for (int it = Lp[i]; it < Lp[i + 1]; it++) {
			int j = Lj[it];
			if (qinv[i] == j) {
				q[j] = i;
				spasm_swap(Lj, it, Lp[i]);
				spasm_swap(Lx, it, Lp[i]);
				break;
			}
		}
		assert(Lj[Lp[i]] == qinv[i]);
		assert(Lx[Lp[i]] == 1);
	}

	/* allocate result and workspace */
	spasm *K = spasm_csr_alloc(n - m, n, L->nzmax, prime, SPASM_WITH_NUMERICAL_VALUES);
	int *xj = malloc(3 * m * sizeof(*xj));
	for (int i = 0; i < 3*m; i++)
		xj[i] = 0;
	spasm_GFp * x = malloc(m * sizeof(*x));

	int nz = 0;
	int *Kp = K->p;
	int *Kj = K->j;
	spasm_GFp *Kx = K->x;
	int Ki = 0;
	for (int i = 0; i < n; i++) {
		if (qinv[i] >= 0)
			continue;   /* pivotal row */
		int top = spasm_sparse_forward_solve(L, L, i, xj, x, q);

		/* enlarge K if necessary */
		if (nz + m - top + 1 > K->nzmax) {
			spasm_csr_realloc(K, 2 * K->nzmax + m - top + 1);
			Kp = K->p;
			Kj = K->j;
			Kx = K->x;
		}
		
		/* finalize previous row of K */
		Kp[Ki++] = nz;

		for (int it = top; it < m; it++) {
			int j = xj[it];
			if (x[j] != 0) {
				Kj[nz] = q[j];
				Kx[nz] = x[j];
				nz++;
			}
		}
		/* add the extra entry for row i */
		Kj[nz] = i;
		Kx[nz] = prime - 1;
		nz++;
	}

	/* finalize last row */
	Kp[Ki] = nz;

	/* cleanup and return */
	free(xj);
	free(x);
	free(q);
	spasm_free_LU(LU);
	spasm_csr_free(L);
	return K;
}
