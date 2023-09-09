#include <stdlib.h>
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
	// const int *Lj = L->j;
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
		spasm_GFp backup = x[i];
		spasm_scatter(L, i, prime - x[i], b);
		x[i] = backup;
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
	int prime = U->prime;
	for (int i = 0; i < n; i++)
		x[i] = 0;

	for (int i = 0; i < n; i++) {
		int j = (q != SPASM_IDENTITY_PERMUTATION) ? q[i] : i;
		if (b[j] != 0) {
			/* check diagonal entry */
			// const spasm_GFp diagonal_entry = Ux[Up[i]];
			// assert(diagonal_entry == 1);

			/* axpy - inplace */
			x[i] = b[j];
			spasm_GFp backup = x[i];
			spasm_scatter(U, i, prime - x[i], b);
			x[i] = backup;
			b[j] = 0;
		}
	}
	for (int i = 0; i < m; i++) 
		if (b[i] != 0)
			return SPASM_NO_SOLUTION;
	return SPASM_SUCCESS;
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
int spasm_sparse_triangular_solve(const spasm *U, const spasm *B, int k, int *xj, spasm_GFp * x, const int *qinv)
{
	int m = U->m;
	int prime = U->prime;
	// const i64 *Bp = B->p;
	// const int *Bj = B->j;
	// const spasm_GFp *Bx = B->x;

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
		int i = (qinv != NULL) ? (qinv[j]) : j;
		if (i < 0)
			continue;

		/* the pivot entry on row i is 1, so we just have to multiply by -x[j] */
		spasm_GFp backup = x[j];
		spasm_scatter(U, i, prime - x[j], x);
		assert(x[j] == 0);
		x[j] = backup;
	}
	return top;
}

/*
 * Solves X * U == B
 * Serious similarity with spasm_schur, spasm_rref, spasm_kernel, ...
 */
spasm *spasm_trsm(const spasm *U, const int *qinv, const spasm *B)
{
	assert(U->m == B->m);
	int m = B->m;
	int n = B->n;
	spasm *X = spasm_csr_alloc(n, m, spasm_nnz(B), B->prime, SPASM_WITH_NUMERICAL_VALUES);
	i64 *Xp = X->p;
	int *Xj = X->j;
	spasm_GFp *Xx = X->x;
	i64 nnz = 0;      /* nnz in X at the moment */
	int Xn = 0;       /* #rows in X at the moment */
	int writing = 0;
	double start = spasm_wtime();

	#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
		int *xj = spasm_malloc(3 * m * sizeof(int));
		for (int j = 0; j < 3*m; j++)
			xj[j] = 0;
		int tid = spasm_get_thread_num();

		#pragma omp for schedule(guided)
		for (int i = 0; i < n; i++) {
			int top = spasm_sparse_triangular_solve(U, B, i, xj, x, qinv);

			int row_nnz = 0;             /* #nz coefficients in the row */
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (x[j] == 0)
					continue;
				assert(qinv[j] >= 0);   /* otherwise, solution does not exist */
				row_nnz += 1;
			}

			int local_i;
			i64 local_nnz;
			#pragma omp critical(schur_complement)
			{
				/* enough room in X? */
				if (nnz + row_nnz > X->nzmax) {
					/* wait until other threads stop writing into it */
					#pragma omp flush(writing)
					while (writing > 0) {
						#pragma omp flush(writing)
					}
					spasm_csr_realloc(X, 2 * X->nzmax + m);
					Xj = X->j;
					Xx = X->x;
				}
				/* save row Xn */
				local_i = Xn;
				Xn += 1;
				local_nnz = nnz;
				nnz += row_nnz;
				#pragma omp atomic update
				writing += 1;    /* register as a writing thread */
			}
			
			/* write the new row in X */
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (x[j] == 0)
					continue;
				Xj[local_nnz] = qinv[j];
				Xx[local_nnz] = x[j];
				local_nnz += 1;
			}
			Xp[local_i + 1] = local_nnz;

			#pragma omp atomic update
			writing -= 1;        /* unregister as a writing thread */

			if (tid == 0) {
				fprintf(stderr, "\r[trsm] %d/%d [%" PRId64 " nz]", Xn, n, nnz);
				fflush(stderr);
			}
		}
		free(x);
		free(xj);
	}
	/* finalize X */
	spasm_csr_realloc(X, -1);
	double density = 1.0 * nnz / (1.0 * m * n);
	fprintf(stderr, "\r[trsm] %d * %d [%" PRId64 " nz / density= %.3f], %.1fs\n", n, m, nnz, density, spasm_wtime() - start);
	return X;
}