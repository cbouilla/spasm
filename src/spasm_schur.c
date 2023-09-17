#include <stdlib.h>
#include <assert.h>
#include <err.h>

#include "spasm.h"

/*
 * Samples R rows at random in the schur complement of (P*A)[0:n] w.r.t. U, and return the average density.
 * qinv locates the pivots in U.
 */
double spasm_schur_estimate_density(const spasm *A, const int *p, int n, const spasm *U, const int *qinv, int R)
{
	if (n == 0)
		return 0;
	int m = A->m;
	i64 nnz = 0;
	if (n == 0)
		return 0;

	#pragma omp parallel
	{
		/* per-thread scratch space */
		spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * m * sizeof(*xj));
		for (int j = 0; j < 3 * m; j++)
			xj[j] = 0;

		#pragma omp for reduction(+:nnz) schedule(dynamic)
		for (int i = 0; i < R; i++) {
			/* pick a random non-pivotal row in A */
			int inew = p[rand() % n];
			int top = spasm_sparse_triangular_solve(U, A, inew, xj, x, qinv);
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if ((qinv[j] < 0) && (x[j] != 0))
					nnz += 1;
			}
		}

		free(x);
		free(xj);
	}
	return ((double) nnz) / (m - U->n) / R;
}

/*
 * Computes the Schur complement of (P*A)[0:n] w.r.t. U
 * The pivots need not be the first entries on the rows.
 * The pivots must be unitary.	
 * This returns a sparse representation of S. 
 *
 * qinv describes the location of pivots in U. qinv[j] == i --> pivot on col j is on row i / -1 if none)
 * note that it is possible to have U == A.
 *
 * if L is not NULL, then the corresponding entries will be added
 * It is understood that row i of A corresponds to row p_in[i] of the original matrix.
 * if p_out is not NULL, then row i of the output corresponds to row p_out[i] of the original matrix.
 *
 * If the estimated density is unknown, set it to -1: it will be evaluated
 */
spasm *spasm_schur(const spasm *A, const int *p, int n, const spasm_lu *fact, 
	double est_density, spasm_triplet *L, const int *p_in, int *p_out)
{
	assert(p != NULL);

	int m = A->m;
	const int *qinv = fact->Uqinv;
	int verbose_step = spasm_max(1, n / 1000);
	if (est_density < 0)
		est_density = spasm_schur_estimate_density(A, p, n, fact->U, qinv, 100);
	long long size = (est_density * n) * m;
	if (size > 2147483648)
		errx(1, "Matrix too large (more than 2^31 entries)");
	i64 prime = spasm_get_prime(A);
	spasm *S = spasm_csr_alloc(n, m, size, prime, SPASM_WITH_NUMERICAL_VALUES);
	i64 *Sp = S->p;
	int *Sj = S->j;
	spasm_ZZp *Sx = S->x;
	i64 snz = 0;                               /* nnz in S at the moment */
	int Sn = 0;                                /* #rows in S at the moment */
	i64 lnz = (L != NULL) ? L->nz : 0;         /* nnz in L at the moment */
	int *Li = (L != NULL) ? L->i : NULL;
	int *Lj = (L != NULL) ? L->j : NULL;
	spasm_ZZp *Lx = (L != NULL) ? L->x : NULL;
	int writing = 0;
	double start = spasm_wtime();

	#pragma omp parallel
	{
		/* scratch space for the triangular solver */
		spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * m * sizeof(*xj));
		for (int j = 0; j < 3 * m; j++)
			xj[j] = 0;
		int tid = spasm_get_thread_num();

		#pragma omp for schedule(dynamic, verbose_step)
		for (int i = 0; i < n; i++) {
			int inew = p[i];
			int top = spasm_sparse_triangular_solve(fact->U, A, inew, xj, x, qinv);

			int row_snz = 0;             /* #nz coefficients in the row of S */
			int row_lnz = 0;             /* #nz coefficients in the row of L */
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (x[j] == 0)
					continue;
				if (qinv[j] < 0)
					row_snz += 1;
				else
					row_lnz += 1;
			}

			int local_i;
			i64 local_snz, local_lnz;
			#pragma omp critical(schur_complement)
			{
				/* enough room in S? */
				if (snz + row_snz > S->nzmax) {
					/* wait until other threads stop writing into it */
					#pragma omp flush(writing)
					while (writing > 0) {
						#pragma omp flush(writing)
					}
					spasm_csr_realloc(S, 2 * S->nzmax + m);
					Sj = S->j;
					Sx = S->x;
				}
				/* save row Sn */
				local_i = Sn;
				Sn += 1;
				local_snz = snz;
				snz += row_snz;

				if (L != NULL && lnz + row_lnz > L->nzmax) {
					/* wait until other threads stop writing into it */
					#pragma omp flush(writing)
					while (writing > 0) {
						#pragma omp flush(writing)
					}
					spasm_triplet_realloc(L, 2 * L->nzmax + m);
					Li = L->i;
					Lj = L->j;
					Lx = L->x;
				}
				local_lnz = lnz;
				lnz += row_lnz;

				#pragma omp atomic update
				writing += 1;    /* register as a writing thread */
			}
			
			/* write the new row in L / S */
			int i_out = (p_in != NULL) ? p_in[inew] : inew;
			if (p_out != NULL)
				p_out[local_i] = i_out;

			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (x[j] == 0)
					continue;
				if (qinv[j] < 0) {
					Sj[local_snz] = j;
					Sx[local_snz] = x[j];
					local_snz += 1;
				} else if (L != NULL) {
					Li[local_lnz] = i_out;
					Lj[local_lnz] = qinv[j];
					Lx[local_snz] = x[j];
					local_lnz += 1;
				}
			}
			Sp[local_i + 1] = local_snz;

			#pragma omp atomic update
			writing -= 1;        /* unregister as a writing thread */

			if (tid == 0 && (i % verbose_step) == 0) {
				double density =  1.0 * snz / (1.0 * m * Sn);
				fprintf(stderr, "\rSchur complement: %d/%d [%" PRId64 " nz / density= %.3f]", Sn, n, snz, density);
				fflush(stderr);
			}
		}
		free(x);
		free(xj);
	}
	/* finalize S */
	spasm_csr_realloc(S, -1);
	double density = 1.0 * snz / (1.0 * m * n);
	fprintf(stderr, "\rSchur complement: %d * %d [%" PRId64 " nz / density= %.3f], %.1fs\n", n, m, snz, density, spasm_wtime() - start);
	return S;
}

static void prepare_q(int m, const int *qinv, int *q)
{
	int i = 0;
	for (int j = 0; j < m; j++)
		if (qinv[j] < 0) {
			q[i] = j;
			i += 1;
		}
}

/*
 * Computes the dense schur complement of (P*A)[0:n] w.r.t. U. 
 * S must be preallocated of dimension n * (A->m - U->n)
 * zero rows are not written to S.
 * return the number of rows actually written to S.
 * S implicitly has dimension k x (m - npiv), row major, lds == m-npiv.
 * q must be preallocated of size at least (m - U->n).
 * on output, q sends columns of S to non-pivotal columns of A
 */
int spasm_schur_dense(const spasm *A, const int *p, int n, const spasm *U, const int *qinv, double *S, int *q)
{
	assert(p != NULL);
	int m = A->m;
	int Sm = m - U->n;                                   /* #columns of S */
	prepare_q(m, qinv, q);
	fprintf(stderr, "[schur/dense] dimension %d x %d...\n", n, Sm);
	double start = spasm_wtime();
	int verbose_step = spasm_max(1, n / 1000);
	i64 r = 0;

	#pragma omp parallel
	{
		/* per-thread scratch space */
		spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * m * sizeof(*xj));
		for (int j = 0; j < 3 * m; j++)
			xj[j] = 0;
		int tid = spasm_get_thread_num();

		#pragma omp for schedule(dynamic, verbose_step)
		for (int k = 0; k < n; k++) {
			int i = p[k];          /* corresponding row of A */
			
			/* eliminate known sparse pivots, put result in x */
			for (int j = 0; j < Sm; j++)
				x[q[j]] = 0;
			int top = spasm_sparse_triangular_solve(U, A, i, xj, x, qinv);
			if (top == m)
				continue;       /* skip empty row */
			
			/* acquire the next row of S */
			i64 t;
			#pragma omp atomic capture
			{ t = r; r += 1; }

			/* gather x into S[t] */
			for (i64 j = 0; j < Sm; j++)
				S[t * Sm + j] = x[q[j]];

			/* verbosity */
			if (tid == 0 && (i % verbose_step) == 0) {
				fprintf(stderr, "\r[schur/dense] %" PRId64 "/%d", r, n);
				fflush(stderr);
			}
		}
		free(x);
		free(xj);
	}
	fprintf(stderr, "\n[schur/dense] finished in %.1fs, rank <= %" PRId64 "\n", spasm_wtime() - start, r);
	return r;
}


/*
 * Computes N random linear combinations rows of the Schur complement of (P*A)[0:n] w.r.t. U.
 * Assumes that pivots are first entries of the row
 * if w > 0, take random linear combinations of subsets of w rows, 
 *    otherwise, take random linear combinations of all the rows
 * S must be preallocated of dimension N * (A->m - U->n)
 * S implicitly has dimension N x (m - npiv), row major, lds == m-npiv.
 * q must be preallocated of size at least (m - U->n).
 * on output, q sends columns of S to non-pivotal columns of A
 */
void spasm_schur_dense_randomized(const spasm *A, const int *p, int n, const spasm *U, const int *qinv, double *S, int *q, int N, int w)
{
	assert(p != NULL);
	assert(n > 0);
	int m = A->m;
	int Sm = m - U->n;
	const i64 *Up = U->p;
	const int *Uj = U->j;
	prepare_q(m, qinv, q);
	fprintf(stderr, "[schur/dense/random] dimension %d x %d, weight %d...\n", N, Sm, w);
	double start = spasm_wtime();
	int verbose_step = spasm_max(1, N / 1000);
	i64 prime = spasm_get_prime(A);

	#pragma omp parallel
	{
		/* per-thread scratch space */
		spasm_ZZp *x = spasm_malloc(m * sizeof(*x));

		#pragma omp for schedule(dynamic, verbose_step)
		for (i64 k = 0; k < N; k++) {
			for (int j = 0; j < m; j++)
				x[j] = 0;
			if (w <= 0) {
				/* x <--- random linear combinations of all rows */
				for (int i = 0; i < n; i++) {
					int inew = p[i];
					int coeff = spasm_ZZp_init(A->field, rand());
					spasm_scatter(A, inew, coeff, x);
				}
			} else {
				for (int i = 0; i < w; i++) {
					int inew = p[rand() % n];
					int coeff = (i == 0) ? 1 : spasm_ZZp_init(A->field, rand());
					spasm_scatter(A, inew, coeff, x);
				}
			}

			/* eliminate known sparse pivots */
			for (int i = 0; i < U->n; i++) {
				int j = Uj[Up[i]];
				if (x[j] == 0)
					continue;
				spasm_scatter(U, i, prime - x[j], x);
			}
			
			/* gather x into S[k] */
			for (i64 j = 0; j < Sm; j++)
				S[k * Sm + j] = x[q[j]];

			/* verbosity */
			if ((k % verbose_step) == 0) {
				fprintf(stderr, "\r[schur/dense/random] %" PRId64 "/%d", k, N);
				fflush(stderr);
			}
		}
		free(x);
	}
	fprintf(stderr, "\n[schur/dense/random] finished in %.1fs\n", spasm_wtime() - start);
}