#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

/*
 * Solving triangular systems, dense RHS
 */


/*
 * Solve x.L = b with dense b and x.
 * x must have size n (#rows of L) and b must have size m (#cols of L)
 * b is destroyed
 * 
 * L is assumed to be (permuted) lower-triangular, with non-zero diagonal.
 * 
 * p[j] == i indicates if the "diagonal" entry on column j is on row i
 * 
 */
void spasm_dense_back_solve(const struct spasm_csr *L, spasm_ZZp *b, spasm_ZZp *x, const int *p)
{
	int n = L->n;
	int r = L->m;
	const i64 *Lp = L->p;
	const int *Lj = L->j;
	const spasm_ZZp *Lx = L->x;
	
	for (int i = 0; i < n; i++)
		x[i] = 0;

	for (int j = r - 1; j >= 0; j--) {
		int i = (p != NULL) ? p[j] : j;
		assert(0 <= i);
		assert(i < n);

		/* scan L[i] to locate the "diagonal" entry on column j */
		spasm_ZZp diagonal_entry = 0;
		for (i64 px = Lp[i]; px < Lp[i + 1]; px++)
			if (Lj[px] == j) {
				diagonal_entry = Lx[px];
				break; 
			}
		assert(diagonal_entry != 0);

		/* axpy - inplace */
		spasm_ZZp alpha = spasm_ZZp_inverse(L->field, diagonal_entry);
		x[i] = spasm_ZZp_mul(L->field, alpha, b[j]);
		spasm_ZZp backup = x[i];
		spasm_scatter(L, i, -x[i], b);
		x[i] = backup;
	}
}

/*
 * Solve x.U = b with dense x, b.
 * 
 * b is destroyed on output
 * 
 * U is (petmuted) upper-triangular with unit diagonal.
 * q[i] == j    means that the pivot on row i is on column j (this is the inverse of the usual qinv). 
 *
 * returns True if a solution was found;
 */
bool spasm_dense_forward_solve(const struct spasm_csr *U, spasm_ZZp *b, spasm_ZZp *x, const int *q)
{
	int n = U->n;
	int m = U->m;
	assert(n <= m);

	for (int i = 0; i < n; i++)
		x[i] = 0;

	for (int i = 0; i < n; i++) {
		int j = (q != NULL) ? q[i] : i;
		
		if (b[j] == 0)
			continue;

		/* eliminate b[j] */
		x[i] = b[j];
		spasm_scatter(U, i, -b[j], b);
		assert(b[j] == 0);
	}
	for (int j = 0; j < m; j++)   /* check that everything has been eliminated */
		if (b[j] != 0)
			return 0;
	return 1;
}

/*
 * solve x * U = B[k], where U is (permuted) triangular (either upper or lower).
 *
 * x must have size m (#columns of U); it does not need to be initialized.
 * xj must be preallocated of size 3*m and zero-initialized (it remains OK)
 * qinv locates the pivots in U.
 *
 * On output, the solution is scattered in x, and its pattern is given in xj[top:m].
 * The precise semantics is as follows. Define:
 *         x_a = { j in [0:m] : qinv[j] < 0 }
 *         x_b = { j in [0:m] : qinv[j] >= 0 }
 * Then x_b * U + x_a == B[k].  It follows that x * U == y has a solution iff x_a is empty.
 * 
 * top is the return value
 *
 * This does not require the pivots to be the first entry of the row.
 * This requires that the pivots in U are all equal to 1. 
 */
int spasm_sparse_triangular_solve(const struct spasm_csr *U, const struct spasm_csr *B, int k, int *xj, spasm_ZZp * x, const int *qinv)
{
	int m = U->m;
	assert(qinv != NULL);
	// const i64 *Bp = B->p;
	// const int *Bj = B->j;
	// const spasm_ZZp *Bx = B->x;

	/* compute non-zero pattern of x --- xj[top:m] = Reach(U, B[k]) */
	int top = spasm_reach(U, B, k, m, xj, qinv);

	/* clear x and scatter B[k] into x*/
	for (int px = top; px < m; px++) {
		int j = xj[px];
		x[j] = 0;
	}
	spasm_scatter(B, k, 1, x);
	// for (i64 px = Bp[k]; px < Bp[k + 1]; px++) {
	// 	int j = Bj[px];
	// 	x[j] = Bx[px];
	// }

	/* iterate over the (precomputed) pattern of x (= the solution) */
	for (int px = top; px < m; px++) {
		int j = xj[px];          /* x[j] is generically nonzero, (i.e., barring numerical cancelation) */

		/* locate corresponding pivot if there is any */
		int i = qinv[j];
		if (i < 0)
			continue;

		/* the pivot entry on row i is 1, so we just have to multiply by -x[j] */
		spasm_ZZp backup = x[j];
		spasm_scatter(U, i, -x[j], x);
		assert(x[j] == 0);
		x[j] = backup;
	}
	return top;
}