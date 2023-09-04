#include <assert.h>

#include "spasm.h"

/*
 * Solving triangular systems, dense RHS
 */


/*
 * dense backwards substitution solver. Solve x . L = b where x and b are
 * dense.
 * 
 * b is undefined on output
 * 
 * L is assumed to be lower-triangular, with non-zero diagonal.
 * 
 * The diagonal entry is the **last** of each row. More precisely, L[j,j] is Lx[
 * Lp[j+1] - 1 ]
 * 
 * p[j] == i indicates if the "diagonal" entry on column j is on row i
 * 
 */
void spasm_dense_back_solve(const spasm *L, spasm_GFp *b, spasm_GFp *x, const int *p)
{
	/* check inputs */
	assert(b != NULL);
	assert(x != NULL);
	assert(L != NULL);
	int n = L->n;
	int m = L->m;
	const i64 *Lp = L->p;
	const int *Lj = L->j;
	const spasm_GFp *Lx = L->x;
	int prime = L->prime;

	for (int i = 0; i < n; i++)
		x[i] = 0;

	for (int j = m - 1; j >= 0; j--) {
		int i = (p != SPASM_IDENTITY_PERMUTATION) ? p[j] : j;

		/* pivot on the j-th column is on the i-th row */
		const spasm_GFp diagonal_entry = Lx[Lp[i + 1] - 1];

		/* axpy - inplace */
		x[i] = (b[j] * spasm_GFp_inverse(diagonal_entry, prime)) % prime;
		spasm_scatter(Lj, Lx, Lp[i], Lp[i + 1] - 1, prime - x[i], b, prime);
	}
}

/*
 * dense forwards substitution solver. Solve x . U = b where x and b are
 * dense.
 * 
 * b is undefined on output
 * 
 * U is upper-triangular
 * 
 * Assumption : the diagonal entry is always present, is always != 0.
 * 
 * The diagonal entry is the first one of each row. More precisely, U[i,i] is
 * Ux[ Up[i] ]
 * 
 * if q != SPASM_IDENTITY_PERMUTATION, then q[i] indicates the column on which
 * the i-th row pivot is.
 * 
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_dense_forward_solve(const spasm * U, spasm_GFp *b, spasm_GFp *x, const int *q)
{
	/* check inputs */
	assert(b != NULL);
	assert(x != NULL);
	assert(U != NULL);
	int n = U->n;
	int m = U->m;
	assert(n <= m);
	const i64 *Up = U->p;
	const int *Uj = U->j;
	const spasm_GFp *Ux = U->x;
	int prime = U->prime;
	for (int i = 0; i < n; i++)
		x[i] = 0;

	for (int i = 0; i < n; i++) {
		int j = (q != SPASM_IDENTITY_PERMUTATION) ? q[i] : i;
		if (b[j] != 0) {
			/* check diagonal entry */
			const spasm_GFp diagonal_entry = Ux[Up[i]];
			assert(diagonal_entry == 1);

			/* axpy - inplace */
			x[i] = b[j];
			spasm_scatter(Uj, Ux, Up[i] + 1, Up[i + 1], prime - x[i], b, prime);
			b[j] = 0;
		}
	}
	for (int i = 0; i < m; i++) 
		if (b[i] != 0)
			return SPASM_NO_SOLUTION;
	return SPASM_SUCCESS;
}



/*************** Triangular solving with sparse RHS
 *
 * solve x * U = B[k], where U is (permuted) upper triangular.
 *
 * x has size m (number of columns of U, paradoxically).
 *
 * when this function returns, the solution is scattered in x, and its pattern
 * is given in xj[top : m].
 *
 * xj must be zero-initialized on the first call (and it stays OK)
 * x does not need to be initialized.
 *
 * top is the return value.
 */
int spasm_sparse_forward_solve(const spasm *U, const spasm *B, int k, int *xj, spasm_GFp * x, const int *qinv)
{
	int m = U->m;
	const i64 *Up = U->p;
	const int *Uj = U->j;
	const spasm_GFp *Ux = U->x;
	int prime = U->prime;
	const i64 *Bp = B->p;
	const int *Bj = B->j;
	const spasm_GFp *Bx = B->x;

	/* compute non-zero pattern of x --- xj[top:m] = Reach(U, B[k]) */
	int top = spasm_reach(U, B, k, m, xj, qinv);

	/* clear x and scatter B[k] into x*/
	for (int px = top; px < m; px++) {
		int j = xj[px];
		x[j] = 0;
	}
	for (i64 px = Bp[k]; px < Bp[k + 1]; px++) {
		int j = Bj[px];
		x[j] = Bx[px];
	}

	/* iterate over the (precomputed) pattern of x (= the solution) */
	for (int px = top; px < m; px++) {
		/* x[j] is nonzero */
		int j = xj[px];

		/* locate corresponding pivot if there is any */
		int i = (qinv != NULL) ? (qinv[j]) : j;
		if (i < 0)
			continue;

		/* the pivot entry on row i is 1, so we just have to multiply by -x[j] */
		assert(Ux[Up[i]] == 1);
		assert(i < U->n);
		spasm_scatter(Uj, Ux, Up[i] + 1, Up[i + 1], prime - x[j], x, prime);
	}
	return top;
}
