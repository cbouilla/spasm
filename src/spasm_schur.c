#include <stdlib.h>
#include <assert.h>
#include <err.h>

#include "spasm.h"

/* make pivotal rows of A unitary */
void spasm_make_pivots_unitary(spasm * A, const int *p, const int npiv)
{
	int prime = A->prime;
	int *Ap = A->p;
	spasm_GFp *Ax = A->x;

	#pragma omp parallel for
	for (int i = 0; i < npiv; i++) {
		int inew = p ? p[i] : i;
		spasm_GFp diag = Ax[Ap[inew]];
		if (diag == 1)
			continue;
		spasm_GFp alpha = spasm_GFp_inverse(diag, prime);
		for (int px = Ap[inew]; px < Ap[inew + 1]; px++)
			Ax[px] = (alpha * Ax[px]) % prime;
	}
}

/* keep only non-pivotal columns and renumber them starting from 0 */
void spasm_stack_nonpivotal_columns(spasm *A, int *qinv)
{
	int n = A->n;
	int m = A->m;
	int *Ap = A->p;
	int *Aj = A->j;
	int *q = spasm_malloc(m * sizeof(int));
	int k = 0;
	for (int j = 0; j < m; j++)
		q[j] = (qinv[j] < 0) ? k++ : -1;

	#pragma omp parallel for
	for (int px = 0; px < Ap[n]; px++) {
		Aj[px] = q[Aj[px]];
		assert(Aj[px] >= 0);
	}
	A->m = k;
}


/*
 * Computes the Schur complement of A w.r.t. U
 * Process rows p[npiv:n] of A.
 * non-pivotal rows are p[npiv:n]
 * The pivots must be the first entries on the rows.
 * The pivots must be unitary.	
 * This returns a sparse representation of S. 
 *
 * the Schur complement is computed from rows p[npiv:n] of A
 * qinv describes the location of pivots in U. qinv[j] == i --> pivot on col j is on row i / -1 if none)
 * note that it is possible to have U == A.
 *
 * If the estimated density is unknown, set it to -1: it will be evaluated
 */
spasm *spasm_schur(const spasm * A, const int *p, int npiv, const spasm *U, const int *qinv, double est_density, int keep_L, int *p_out)
{
	assert(!keep_L); /* option presently unsupported */
	
	const int n = A->n;
	const int m = A->m;
	const int Sm = m;
	const int Sn = n - npiv;
	const int verbose_step = spasm_max(1, n / 1000);
	
	if (est_density < 0)
		est_density = spasm_schur_probe_density(A, p, npiv, U, qinv, 100);
	long long size = (est_density * Sn) * Sm;
	if (size > 2147483648)
		errx(1, "Matrix too large (more than 2^31 entries)");
	spasm *S = spasm_csr_alloc(Sn, Sm, (est_density*Sn)*Sm, A->prime, SPASM_WITH_NUMERICAL_VALUES);
	int *Sp = S->p;
	int *Sj = S->j;
	spasm_GFp *Sx = S->x;
	int snz = 0;
	int k = 0;
	int writing = 0;
	double start = spasm_wtime();

	#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
		int *xj = spasm_malloc(3 * m * sizeof(int));
		spasm_vector_zero(xj, 3 * m);
		int tid = spasm_get_thread_num();
		int row_snz, row_k, row_px;

		#pragma omp for schedule(dynamic, verbose_step)
		for (int i = npiv; i < n; i++) {
			int inew = p[i];
			int top = spasm_sparse_forward_solve(U, A, inew, xj, x, qinv);

			/* count surviving coefficients */
			row_snz = 0;
			for (int px = top; px < m; px++) {
				const int j = xj[px];
				if ((keep_L || (qinv[j] < 0)) && (x[j] != 0))
					row_snz++;
			}

			#pragma omp critical(schur_complement)
			{
				/* enough room in S? */
				if (snz + row_snz > S->nzmax) {
					/* wait until other threads stop writing into it */
					#pragma omp flush(writing)
					while (writing > 0) {
						#pragma omp flush(writing)
					}
					spasm_csr_realloc(S, 2 * S->nzmax + Sm);
					Sj = S->j;
					Sx = S->x;
				}
				/* save row k */
				row_k = k++;
				row_px = snz;
				snz += row_snz;
				#pragma omp atomic update
				writing++;
			}
			if (p_out)
				p_out[row_k] = inew;
			
			/* write the new row in S */
			Sp[row_k] = row_px;
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if ((keep_L || (qinv[j] < 0)) && (x[j] != 0)) {
					Sj[row_px] = j;
					Sx[row_px++] = x[j];
				}
			}

			#pragma omp atomic update
			writing--;

			if (tid == 0 && (i % verbose_step) == 0) {
				double density =  1.0 * snz / (1.0 * Sm * k);
				fprintf(stderr, "\rSchur complement: %d/%d [%d NNZ / density= %.3f]", k, Sn, snz, density);
				fflush(stderr);
			}
		}
		free(x);
		free(xj);
	}
	/* finalize S */
	Sp[Sn] = snz;
	spasm_csr_realloc(S, -1);
	double density = 1.0 * snz / (1.0 * Sm * Sn);
	fprintf(stderr, "\rSchur complement: %d * %d [%d NNZ / density= %.3f], %.1fs\n", Sn, Sm, snz, density, spasm_wtime() - start);
	return S;
}


/** Samples R rows at random in the schur complement of A w.r.t. the pivots in p[0:n_pivots],
* and return the number that are non-zero (these rows of A are linearly independent from the pivots).
* The pivots must be unitary.
*/
double spasm_schur_probe_density(const spasm * A, const int *p, int npiv, const spasm *U, const int *qinv, int R)
{
	int nnz = 0;
	int m = A->m;
	int n = A->n;
	if (m == npiv || n == npiv)
		return 0.0;

	#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
		int *xj = spasm_malloc(3 * m * sizeof(int));
		spasm_vector_zero(xj, 3 * m);

		#pragma omp for reduction(+:nnz) schedule(dynamic)
		for (int i = 0; i < R; i++) {
			/* pick a random row in S, check if non-zero */
			int inew = p[npiv + (rand() % (n - npiv))];
			int top = spasm_sparse_forward_solve(U, A, inew, xj, x, qinv);
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (qinv[j] < 0 && x[j] != 0)
					nnz++;
			}
		}
		free(x);
		free(xj);
	}
	return ((double) nnz) / (m - npiv) / R;
}

/*
 * computes the rank of the schur complement, but not the schur complement
 * itself. The pivots must be unitary.
 */
int spasm_schur_rank_old(spasm * A, const int *p, const int *qinv, const int npiv) {
	int Sm, m, n, k, r, prime, step, threads, searched, prev_r;
	int *q, *Ap, *Aj;
	double start;
	spasm_GFp *Ax;

	n = A->n;
	m = A->m;
	Ap = A->p;
	Aj = A->j;
	Ax = A->x;
	prime = A->prime;

	/* Get Workspace */
	Sm = m - npiv;
	q = spasm_malloc(Sm * sizeof(int));

	/* q sends columns of S to non-pivotal columns of A */
	k = 0;
	for (int j = 0; j < m; j++)
		if (qinv[j] < 0)
			q[k++] = j;

	spasm_dense_lu *U = spasm_dense_LU_alloc(Sm, A->prime);

	/* ---- compute Schur complement ----- */
	fprintf(stderr, "rank of dense schur complement...\n");

	start = spasm_wtime();
	r = 0;
	step = 1;
	k = 0;
	searched = 0;
	prev_r = 0;
	threads = spasm_get_num_threads();

	#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
		spasm_GFp *y = spasm_malloc(Sm * sizeof(spasm_GFp));
		int gain;

		while (step <= (1 << 16)) {	/* <--- tweak-me */
			double it_start = spasm_wtime();
			prev_r = r;

			/* random linear combination */
			spasm_vector_zero(x, m);
			for (int i = 0; i < step; i++) {
				int inew = p[npiv + (rand() % (n - npiv))];
				spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], 1 + (rand() % (prime - 1)), x, prime);
			}
			spasm_eliminate_sparse_pivots(A, npiv, p, x);
			for (int j = 0; j < Sm; j++)	/* gather into y */
				y[j] = x[q[j]];

			#pragma omp atomic update
			r += spasm_dense_LU_process(U, y);

			/* this is a barrier */
			#pragma omp single
			{
				fprintf(stderr, "\rSchur rank: %d [%.1fs] -- current rank = %d / step = %d", k, spasm_wtime() - it_start, r, step);
				fflush(stderr);

				k++;
				searched += threads * step;
				gain = r - prev_r;

				if (gain < threads)
					step *= 2;
				else
					step = spasm_max(1, step / 2);
			}
		}

		#pragma omp single
		{
			int final_bad = 0;
			k = 0;
			fprintf(stderr, "\n");

			while (final_bad < 3) {
				double it_start = spasm_wtime();
				for (int i = npiv; i < n; i++) {
					int inew = p[i];
					spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, x, prime);
				}
				spasm_eliminate_sparse_pivots(A, npiv, p, x);
				for (int j = 0; j < Sm; j++)
					y[j] = x[q[j]];
				int new = spasm_dense_LU_process(U, y);
				r += new;
				final_bad += 1 - new;
				k++;
				fprintf(stderr, "\rSchur rank: %d [%.1fs] -- current rank = %d / final", k, spasm_wtime() - it_start, r);
				fflush(stderr);
			}
		}
		free(x);
		free(y);
	}
	fprintf(stderr, "\n[schur/rank] Time: %.1fs\n", spasm_wtime() - start);

	free(q);
	spasm_dense_LU_free(U);
	return r;
}


/*
 * Computes the dense schur complement. Generates rows [lo:hi] of the schur complement. 
 * Usual conditions on pivots apply.
 * returns the number of rows written to S.
 * S implicitly has dimension (hi - lo) x (m - npiv), row major, lds == m-npiv.
 */
int spasm_schur_dense(const spasm *A, const int *p, const int *qinv, const int npiv, 
	              int lo, int hi, double *S)
{
	int n = A->n;
	int m = A->m;
	int Sm = m - npiv;                                   /* #columns of S */
	int Sn = hi - lo;
	assert(0 <= lo);                                     /* validate args */
	assert(lo <= hi);
	assert(hi <= n - npiv);

	/* q sends columns of S to non-pivotal columns of A */
	int *q = spasm_malloc(Sm * sizeof(*q));
	int k = 0;
	for (int j = 0; j < m; j++)
		if (qinv[j] < 0)
			q[k++] = j;

	fprintf(stderr, "Dense schur complement (rows [%d:%d]...\n", lo, hi);
	double start = spasm_wtime();
	int verbose_step = spasm_max(1, Sn / 1000);
	int r = 0;

	#pragma omp parallel
	{
		/* per-thread scratch space */
		spasm_GFp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * m * sizeof(*xj));
		spasm_vector_zero(xj, 3 * m);
		int tid = spasm_get_thread_num();

		#pragma omp for schedule(dynamic, verbose_step)
		for (int k = lo; k < hi; k++) {
			int i = p[npiv + k];   /* corresponding row of A */
			
			/* eliminate known sparse pivots, put result in x */
			for (int j = 0; j < Sm; j++)
				x[q[j]] = 0;
			int top = spasm_sparse_forward_solve(A, A, i, xj, x, qinv);
			if (top == m)
				continue;       /* skip empty row */
			
			/* acquire the next row of S */
			int t;
			#pragma omp atomic capture
			{ t = r; r += 1; }

			/* gather x into S[t] */
			for (int j = 0; j < Sm; j++)
				S[t * Sm + j] = x[q[j]];

			/* verbosity */
			if (tid == 0 && (i % verbose_step) == 0) {
				fprintf(stderr, "\rSchur complement (dense): %d/%d", k, Sn);
				fflush(stderr);
			}
		}
		free(x);
		free(xj);
	}
	fprintf(stderr, "\n[schur/dense] Time: %.1fs\n", spasm_wtime() - start);
	free(q);
	return r;
}

int spasm_schur_rank(const spasm * A, const int *p, const int *qinv, const int npiv)
{
	int n = A->n;
	int m = A->m;
	/* Get Workspace */
	size_t Sm = m - npiv;
	size_t Sn = spasm_min(n - npiv, m - npiv);   /* upper-bound on the rank of the schur complement */
	double *S = spasm_malloc(Sn * Sm * sizeof(*S));

	size_t *Q = spasm_malloc(Sm * sizeof(*Q));

	for (int lo = 0; lo < n - npiv; lo += Sn) {
		int hi = spasm_min(n - npiv, lo + Sn);
		int r = spasm_schur_dense(A, p, qinv, npiv, lo, hi, S);
		fprintf(stderr, "processed A[%d:%d], rank <= %d (out of %zd rows)\n", lo, hi, r, Sn);
		int rr = spasm_dense_echelonize(A->prime, r, Sm, S, m, Q);
		fprintf(stderr, "processed A[%d:%d], actual rank = %d\n", lo, hi, rr);
	}
	return 42;
}