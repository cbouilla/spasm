#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

/*
 * (dense) vector * (sparse) Matrix
 * y <--- x*A + y
 */
void spasm_xApy(const spasm_GFp *x, const spasm *A, spasm_GFp *y)
{
	int n = A->n;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_GFp *Ax = A->x;
	int prime = A->prime;
	for (int i = 0; i < n; i++)
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			y[j] = (y[j] + x[i] * Ax[px]) % prime;
		}
}

/*
 * (sparse) Matrix * (dense) vector
 * y <--- A*x + y
 */
void spasm_Axpy(const spasm *A, const spasm_GFp *x, spasm_GFp *y)
{
	int n = A->n;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_GFp *Ax = A->x;
	int prime = A->prime;
	for (int i = 0; i < n; i++)
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			y[i] = (y[i] + Ax[px] * x[j]) % prime;
		}
}

#if 0
/*
 * (sparse vector) * (sparse) Matrix Compute y = x * M, where x and M are
 * sparse.
 * 
 * The result is scattered in y, its pattern is given by yi. The return value nz
 * is the number of non-zero entries in y.
 */
int spasm_sparse_vector_matrix_prod(const spasm * M, const spasm_GFp * x, const int *xi, int xnz, spasm_GFp * y, int *yi)
{
	/* check inputs */
	i64 Mnz = spasm_nnz(M);
	assert(x != NULL);
	assert(Mnz != 0);

	int m = M->m;
	i64 *Mp = M->p;
	int *Mj = M->j;
	spasm_GFp *Mx = M->x;
	int prime = M->prime;

	/* get workspace, initialize w */
	int *w = spasm_calloc(m, sizeof(*w));

	/* find pattern of result */
	int nz = 0;
	for (i64 k = 0; k < xnz; k++) {
		int i = xi[k];
		for (i64 p = Mp[i]; p < Mp[i + 1]; p++) {
			int j = Mj[p];
			if (w[j] == 0) {
				w[j] = 1;
				yi[nz] = j;
				nz += 1;
			}
		}
	}

	/* form result */
	for (int k = 0; k < xnz; k++) {
		int i = xi[k];
		spasm_scatter(Mj, Mx, Mp[i], Mp[i + 1], x[i], y, prime);
	}

	/* free workspace */
	free(w);
	return nz;
}
#endif