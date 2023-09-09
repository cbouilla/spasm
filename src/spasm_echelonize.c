#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "spasm.h"

/*
 * There are 3 possible finalization strategies
 *   1. GPLU (without computing the sparse schur complement --- GPLU will do it on the fly)
 *   2. dense schur generation + dense factorization
 *   3. low-rank dense (including "tall and skinny") --- pick random linear combinations
 *
 * all echelonization function should have prototype: 
 *     echelonize_XXX(spasm *A, const int *p, int n, spasm *U, int *Uqinv, struct echelonize_opts *opts);
 *
 * passing NULL as the options leads to default choices. 
 * U must be preallocated with sufficiently many rows (e.g. with A->n rows to be ready for the worst case).
 *
 * presently, this does not provide a way to store L; we will do something about this later on.
 */


/* provide sensible defaults */
void spasm_echelonize_init_opts(struct echelonize_opts *opts)
{
	opts->enable_greedy_pivot_search = 1;
	
	opts->enable_tall_and_skinny = 1;
	opts->enable_dense = 1;
	opts->enable_GPLU = 1;

	// options of the main procedure
	opts->min_pivot_proportion = 0.05;
	opts->max_round = 3;
	opts->sparsity_threshold = 0.05;
	opts->tall_and_skinny_ratio = 5;

	opts->dense_block_size = 1000;
	opts->low_rank_ratio = 0.5;
	opts->low_rank_start_weight = -1;
}

bool spasm_echelonize_test_completion(spasm *A, const int *p, int n, spasm *U, int *Uqinv)
{
	int m = A->m;
	int Sm = m - U->n;
	int Sn = ceil(128 / log2(A->prime));
	double *S = spasm_malloc(Sn * Sm * sizeof(*S));
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sp = spasm_malloc(Sm * sizeof(*Sp));       /* for FFPACK */
	fprintf(stderr, "[echelonize/completion] Testing completion with %d random linear combinations (rank %d)\n", Sn, U->n);
	fflush(stderr);
	spasm_schur_dense_randomized(A, p, n, U, Uqinv, S, q, Sn, 0);
	int rr = spasm_ffpack_echelonize(A->prime, Sn, Sm, S, Sm, Sp);
	free(S);
	free(Sp);
	free(q);
	return (rr == 0);
}


/* not dry w.r.t. spasm_LU() */
void spasm_echelonize_GPLU(spasm *A, const int *p, int n, spasm *U, int *qinv, struct echelonize_opts *opts)
{
	int m = A->m;
	int prime = A->prime;
	int verbose_step = spasm_max(1, n / 1000);

	fprintf(stderr, "[echelonize/GPLU] processing matrix of dimension %d x %d\n", n, m);

	/* workspace for triangular solver */
	int *x = spasm_malloc(m * sizeof(spasm_GFp));
	int *xj = spasm_malloc(3 * m * sizeof(int));
	spasm_vector_zero(xj, 3 * m);
	
	/* allocate result */
	i64 *Up = U->p;
	i64 unz = spasm_nnz(U);
	int r = spasm_min(A->n, m);  /* upper-bound on rank */

	/* initialize early abort */
	int rows_since_last_pivot = 0;
	bool early_abort_done = 0;

	/* Main loop : compute L[i] and U[i] */
	int i;
	for (i = 0; i < n; i++) {
		/* test for early abort */
		if (U->n == r) {
			fprintf(stderr, "\n[echelonize/GPLU] full rank reached\n");
			break;
		}
		/* TODO: make these hard-coded values options */
		if (!early_abort_done && rows_since_last_pivot > 10 && (rows_since_last_pivot > (n / 100))) {
			fprintf(stderr, "\n[echelonize/GPLU] testing for early abort...");
			if (spasm_echelonize_test_completion(A, p, n, U, qinv))
				break;
			early_abort_done = 1;
		}

		/* Triangular solve: x * U = A[i] */
		int inew = (p != NULL) ? p[i] : i;
		int top = spasm_sparse_triangular_solve(U, A, inew, xj, x, qinv);

		/* Find pivot column; current poor strategy= choose leftmost */
		int jpiv = -1;	          /* column index of best pivot so far. */
		for (int px = top; px < m; px++) {
			int j = xj[px];        /* x[j] is (generically) nonzero */
			if (x[j] == 0)
				continue;
			if (qinv[j] < 0) {
				/* have found the pivot on row i yet ? */
				if (jpiv == -1 || j < jpiv)
					jpiv = j;
			}
		}
		if (jpiv < 0) {
			/* no pivot found */
			rows_since_last_pivot += 1;
			continue;
		}

		/* ensure enough room in U for an extra row */
		if (unz + m > U->nzmax)
			spasm_csr_realloc(U, 2 * U->nzmax + m);
		int *Uj = U->j;
		spasm_GFp *Ux = U->x;

		/* store new pivotal row into U */
		qinv[jpiv] = U->n;
		i64 old_unz = unz;
		Uj[unz] = jpiv;
		Ux[unz] = 1;
		unz += 1;
		spasm_GFp beta = spasm_GFp_inverse(x[jpiv], prime);
		for (int px = top; px < m; px++) {
			int j = xj[px];
			if (qinv[j] < 0) {
				Uj[unz] = j;
				Ux[unz] = (x[j] * beta) % prime;
				unz += 1;
			}
		}
		U->n += 1;
		Up[U->n] = unz;

		/* TODO: switch to dense */

		/* reset early abort */
		rows_since_last_pivot = 0;
		early_abort_done = 0;
		
		if ((i % verbose_step) == 0) {
			fprintf(stderr, "\r[echelonize/GPLU] %d / %d [|U| = %" PRId64 "] -- current density= (%.3f vs %.3f) --- rank >= %d", i, n, unz, 1.0 * (m - top) / (m), 1.0 * (unz - old_unz) / m, U->n);
			fflush(stderr);
		}
	}
	/* cleanup */
	fprintf(stderr, "\n");
	free(x);
	free(xj);
}


static void dense_update_U(spasm *U, int rr, int Sm, const double *S, const size_t *Sqinv, const int *q, int *Uqinv)
{
	i64 extra_nnz = ((i64) (1 + Sm - rr)) * rr;     /* maximum size increase */
        i64 unz = spasm_nnz(U);
        fprintf(stderr, "[dense update] enlarging U from %" PRId64 " to %" PRId64 " entries\n", unz, unz + extra_nnz);
	spasm_csr_realloc(U, unz + extra_nnz);
	i64 *Up = U->p;
	int *Uj = U->j;
	spasm_GFp *Ux = U->x;
        int prime = U->prime;
        for (i64 i = 0; i < rr; i++) {
                int j = Sqinv[i];   /* column (of S) with the pivot on row i of S; the pivot is implicitly 1 */
        	Uj[unz] = q[j];  /* column of A with the pivot */
        	Ux[unz] = 1;
        	unz += 1;
        	Uqinv[q[j]] = U->n;
        	for (i64 k = rr; k < Sm; k++) {
        		i64 j = Sqinv[k];
        		if (S[i * Sm + k] == 0)
        			continue;   /* don't store zero */
        		Uj[unz] = q[j];
        		Ux[unz] = (prime + (i64) S[i * Sm + k]) % prime;
        		unz += 1;
        	}
        	U->n += 1;
        	Up[U->n] = unz;
        }
        assert(unz == spasm_nnz(U));
}

void spasm_echelonize_dense_lowrank(spasm *A, const int *p, int n, spasm *U, int *Uqinv, struct echelonize_opts *opts)
{
	assert(opts->dense_block_size > 0);
	int m = A->m;
	int Sm = m - U->n;

	i64 size_S = (i64) opts->dense_block_size * (i64) Sm * sizeof(double);
	double *S = spasm_malloc(size_S);
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sp = spasm_malloc(Sm * sizeof(*Sp));       /* for FFPACK */
	double start = spasm_wtime();
	int old_un = U->n;
	int round = 0;
	fprintf(stderr, "[echelonize/dense/low-rank] processing dense schur complement of dimension %d x %d; block size=%d\n", 
		n, Sm, opts->dense_block_size);
	
	/* 
	 * stupid algorithm to decide a starting weight:
	 * - estimate the number of passes as rank_ub / opts->dense_block_size
         * - thus rank_ub * w rows are selected in total
         * - a single row is never selected with proba (1 - 1/n) ** (rank_ub * w)
	 * - choose w such that this is less than 0.01
	 * - this leads to w >= log 0.01 / log (1 - 1/n) / rank_ub
	 * - this is approximately w >= log 0.01 * n / rank_ub
	 */
	int rank_ub = spasm_min(n, Sm);
	int w = (opts->low_rank_start_weight < 0) ? ceil(-log(0.01) * n / rank_ub) : opts->low_rank_start_weight;

	for (;;) {
		/* compute a chunk of the schur complement, then echelonize with FFPACK */
		rank_ub = spasm_min(A->n - U->n, A->m - U->n);
		int Sn = spasm_min(rank_ub, opts->dense_block_size);
		if (Sn <= 0)
			break;		
		fprintf(stderr, "[echelonize/dense/low-rank] Round %d. Weight %d. Processing chunk (%d x %d), |U| = %"PRId64"\n", 
			round, w, Sn, Sm, spasm_nnz(U));
		spasm_schur_dense_randomized(A, p, n, U, Uqinv, S, q, Sn, w);
		int rr = spasm_ffpack_echelonize(A->prime, Sn, Sm, S, Sm, Sp);

		if (rr == 0) {
			if (spasm_echelonize_test_completion(A, p, n, U, Uqinv))
				break;
			fprintf(stderr, "[echelonize/dense/low-rank] Failed termination test; switching to full linear combinations\n");
			w = 0;
			Sn = omp_get_max_threads();
		}
		if (rr < 0.9 * Sn) {
			w *= 2;
			fprintf(stderr, "[echelonize/dense/low-rank] Not enough pivots, increasing weight to %d\n", w);
		}
		dense_update_U(U, rr, Sm, S, Sp, q, Uqinv);
		n -= rr;
		Sm -= rr;
        	round += 1;
        	fprintf(stderr, "[echelonize/dense/low-rank] found %d new pivots (%d new since beginning)\n", rr,  U->n - old_un);
        }
        fprintf(stderr, "[echelonize/dense/low-rank] completed in %.1fs. %d new pivots found\n", spasm_wtime() - start, U->n - old_un);
        free(S);
        free(q);
        free(Sp);
}

/* 
 * the schur complement (on non-pivotal rows of A) w.r.t. U is dense.
 * process (P*A)[0:n]
 */
void spasm_echelonize_dense(spasm *A, const int *p, int n, spasm *U, int *Uqinv, struct echelonize_opts *opts)
{
	assert(opts->dense_block_size > 0);
	int m = A->m;
	int Sm = m - U->n;

	double *S = spasm_malloc(opts->dense_block_size * Sm * sizeof(*S));
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sp = spasm_malloc(Sm * sizeof(*Sp));       /* for FFPACK */
	int processed = 0;
	double start = spasm_wtime();
	int old_un = U->n;
	int round = 0;
	fprintf(stderr, "[echelonize/dense] processing dense schur complement of dimension %d x %d; block size=%d\n", 
		n, Sm, opts->dense_block_size);
	bool lowrank_mode = 0;

	for (;;) {
		/* compute a chunk of the schur complement, then echelonize with FFPACK */
		int rank_ub = spasm_min(A->n - U->n, A->m - U->n);
		int Sn = spasm_min(n - processed, spasm_min(rank_ub, opts->dense_block_size));
		if (Sn <= 0)
			break;
		
		fprintf(stderr, "[echelonize/dense] Round %d. processing S[%d:%d] (%d x %d)\n", round, processed, processed + Sn, Sn, Sm);	
		int r = spasm_schur_dense(A, p, Sn, U, Uqinv, S, q);
		int rr = spasm_ffpack_echelonize(A->prime, r, Sm, S, Sm, Sp);
		
		/* update U */
		dense_update_U(U, rr, Sm, S, Sp, q, Uqinv);

        	/* move on to the next chunk */
        	round += 1;
        	processed += Sn;
        	p += Sn;
        	n -= Sn;
        	Sm = m - U->n;
        	fprintf(stderr, "[echelonize/dense] found %d new pivots\n", rr);

		/* switch to low-rank mode ? */
		if (opts->enable_tall_and_skinny && (rr < opts->low_rank_ratio * Sn)) {
			fprintf(stderr, "[echelonize/dense] Too few pivots; switching to low-rank mode\n");
			lowrank_mode = 1;
			break;
		}
        }
        free(S);
        free(q);
        free(Sp);
        if (lowrank_mode) {
		spasm_echelonize_dense_lowrank(A, p, n, U, Uqinv, opts);
        } else {
        	fprintf(stderr, "[echelonize/dense] completed in %.1fs. %d new pivots found\n", spasm_wtime() - start, U->n - old_un);
        }
}


/*
 * (main entry point)
 * Returns the row echelon form of A. 
 * Initializes Uqinv (must be preallocated of size m [==#columns of A]).
 * Modifies A (permutes entries in rows)
 */
spasm* spasm_echelonize(spasm *A, int *Uqinv, struct echelonize_opts *opts)
{
	struct echelonize_opts default_opts;
	if (opts == NULL) {
		fprintf(stderr, "[echelonize] using default settings\n");
		opts = &default_opts;
		spasm_echelonize_init_opts(opts);
	}
	int n = A->n;
	int m = A->m;
	int prime = A->prime;
	
	int *qinv = spasm_malloc(m * sizeof(int));  /* for pivot search */
	int *p = spasm_malloc(n * sizeof(int));

	spasm *U = spasm_csr_alloc(n, m, spasm_nnz(A), prime, SPASM_WITH_NUMERICAL_VALUES);
	U->n = 0;
	i64 Unz = 0;      /* #entries in U */
	for (int j = 0; j < m; j++)
		Uqinv[j] = -1;

	fprintf(stderr, "[echelonize] Start on %d x %d matrix with %" PRId64 " nnz\n", n, m, spasm_nnz(A));
	double start = spasm_wtime();
	double density = -1;
	int npiv = 0;

	for (int round = 0; round < opts->max_round; round++) {
		fprintf(stderr, "[echelonize] round %d\n", round);
	 	i64 *Ap = A->p;
		int *Aj = A->j;
		spasm_GFp *Ax = A->x;

		/* find structural pivots in A */
		npiv = spasm_find_pivots(A, p, qinv, opts);

		/* compute total pivot nnz and reallocate U if necessary */
		int pivot_nnz = 0;
		for (int k = 0; k < npiv; k++) {
			int i = p[k];
			pivot_nnz += spasm_row_weight(A, i);
		}
		if (spasm_nnz(U) + pivot_nnz > U->nzmax)
			spasm_csr_realloc(U, spasm_nnz(U) + pivot_nnz);

		/* copy pivotal rows to U and make them unitary; update Uqinv */
		i64 *Up = U->p;
		int *Uj = U->j;
		spasm_GFp *Ux = U->x;
		for (int k = 0; k < npiv; k++) {
			int i = p[k];
			int j = Aj[Ap[i]];
			assert(Uqinv[j] < 0);
			assert(qinv[j] == i);
			Uqinv[j] = U->n;
			for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
				Uj[Unz] = Aj[px];
				Ux[Unz] = Ax[px];
				Unz += 1;
			}
			U->n += 1;
			Up[U->n] = Unz;
		}
		assert(Unz <= U->nzmax);
		spasm_make_pivots_unitary(U, SPASM_IDENTITY_PERMUTATION, U->n); /* reprocess previous pivots; shame */

		/* decide whether to move on to the next iteration */
		if (npiv == spasm_min(n, m)) {
			fprintf(stderr, "[echelonize] full rank reached\n");			
			break;
		}

		if (npiv < opts->min_pivot_proportion * spasm_min(n, m)) {
			fprintf(stderr, "[echelonize] not enough pivots found; stopping\n");
			break;     /* not enough pivots found */
		}		
		// if (density > opts->sparsity_threshold && aspect_ratio > opts->tall_and_skinny_ratio) {
		// 	fprintf(stderr, "Schur complement is dense,tall and skinny (#rows / #cols = %.1f)\n", aspect_ratio);
		// 	break;
		// }
		density = spasm_schur_estimate_density(A, p + npiv, n - npiv, U, Uqinv, 100);
		if (density > opts->sparsity_threshold) {
			fprintf(stderr, "[echelonize] Schur complement is dense (estimated %.2f%%)\n", 100 * density);
			break;
		}

		/* compute the next schur complement */
		i64 nnz = (density * (n - npiv)) * (m - npiv);
		char tmp[8];
		spasm_human_format(sizeof(int) * (n - npiv + nnz) + sizeof(spasm_GFp) * nnz, tmp);
		fprintf(stderr, "Schur complement is %d x %d, estimated density : %.2f (%s byte)\n", n - npiv, m - npiv, density, tmp);
		spasm *S = spasm_schur(A, p + npiv, n - npiv, U, Uqinv, density, SPASM_DISCARD_L, NULL);
		if (round > 0)
			spasm_csr_free(A);       /* only if it is not the input argument */
		A = S;
		n = A->n;
		round += 1;
	}
	free(qinv);

	/* finish */
	if (!opts->enable_tall_and_skinny)
		fprintf(stderr, "[echelonize] dense low-rank mode disabled\n");
	if (!opts->enable_dense)
		fprintf(stderr, "[echelonize] regular dense mode disabled\n");
	if (!opts->enable_GPLU)
		fprintf(stderr, "[echelonize] GPLU mode disabled\n");
	
	double aspect_ratio = (double) (n - npiv) / (m - npiv);
	fprintf(stderr, "[echelonize] finishing; density = %.3f; aspect ratio = %.1f\n", density, aspect_ratio);
	if (opts->enable_tall_and_skinny && aspect_ratio > opts->tall_and_skinny_ratio)
		spasm_echelonize_dense_lowrank(A, p + npiv, n - npiv, U, Uqinv, opts);
	else if (opts->enable_dense && density > opts->sparsity_threshold)
		spasm_echelonize_dense(A, p + npiv, n - npiv, U, Uqinv, opts);
	else if (opts->enable_GPLU)
		spasm_echelonize_GPLU(A, p + npiv, n - npiv, U, Uqinv, opts);
	else
		fprintf(stderr, "[echelonize] Cannot finish (no valid method enabled). Incomplete echelonization returned\n");

	free(p);
	fprintf(stderr, "[echelonize] Done in %.1fs. Rank %d, %" PRId64 " nz in basis\n", spasm_wtime() - start, U->n, spasm_nnz(U));
	spasm_csr_resize(U, U->n, m);
	spasm_csr_realloc(U, -1);
	return U;
}