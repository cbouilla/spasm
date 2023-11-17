#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "spasm.h"

/* provide sensible defaults */
void spasm_echelonize_init_opts(struct echelonize_opts *opts)
{
	opts->enable_greedy_pivot_search = 1;
	
	opts->enable_tall_and_skinny = 1;
	opts->enable_dense = 1;
	opts->enable_GPLU = 1;

	// options of the main procedure
	opts->L = 0;
	opts->complete = 0;
	opts->min_pivot_proportion = 0.1;
	opts->max_round = 3;
	opts->sparsity_threshold = 0.05;
	opts->tall_and_skinny_ratio = 5;

	opts->dense_block_size = 1000;
	opts->low_rank_ratio = 0.5;
	opts->low_rank_start_weight = -1;
}

bool spasm_echelonize_test_completion(const struct spasm_csr *A, const int *p, int n, struct spasm_csr *U, int *Uqinv)
{
	/* deal with the easy cases first */
	if (n == 0 || spasm_nnz(A) == 0)
		return 1;
	int m = A->m;
	i64 Sm = m - U->n;
	i64 prime = spasm_get_prime(A);
	spasm_datatype datatype = spasm_datatype_choose(prime);
	i64 Sn = ceil(128 / log2(prime));
	void *S = spasm_malloc(Sn * Sm * spasm_datatype_size(datatype));
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sp = spasm_malloc(Sm * sizeof(*Sp));       /* for FFPACK */
	fprintf(stderr, "[echelonize/completion] Testing completion with %" PRId64" random linear combinations (rank %d)\n", Sn, U->n);
	fflush(stderr);
	spasm_schur_dense_randomized(A, p, n, U, Uqinv, S, datatype, q, Sn, 0);
	int rr = spasm_ffpack_rref(prime, Sn, Sm, S, Sm, datatype, Sp);
	free(S);
	free(Sp);
	free(q);
	return (rr == 0);
}


static void echelonize_GPLU(const struct spasm_csr *A, const int *p, int n, const int *p_in, struct spasm_lu *fact, struct echelonize_opts *opts)
{
	(void) opts;
	assert(p != NULL);
	int m = A->m;
	int r = spasm_min(A->n, m);  /* upper-bound on rank */
	int verbose_step = spasm_max(1, n / 1000);
	fprintf(stderr, "[echelonize/GPLU] processing matrix of dimension %d x %d\n", n, m);
	
	struct spasm_csr *U = fact->U;
	struct spasm_triplet *L = fact->Ltmp;
	int *Uqinv = fact->qinv;
	i64 *Up = U->p;
	i64 unz = spasm_nnz(U);
	i64 lnz = (L != NULL) ? L->nz : 0;
	int *Lp = fact->p;

	/* initialize early abort */
	int rows_since_last_pivot = 0;
	bool early_abort_done = 0;

	/* workspace for triangular solver */
	spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
	int *xj = spasm_malloc(3 * m * sizeof(*xj));
	for (int j = 0; j < 3*m; j++)
		xj[j] = 0;

	/* Main loop : compute L[i] and U[i] */
	int i;
	for (i = 0; i < n; i++) {
		/* test for early abort (if L not needed) */
		if (L == NULL && U->n == r) {
			fprintf(stderr, "\n[echelonize/GPLU] full rank reached\n");
			break;
		}
		/* TODO: make these hard-coded values options */
		if (L == NULL && !early_abort_done && rows_since_last_pivot > 10 && (rows_since_last_pivot > (n / 100))) {
			fprintf(stderr, "\n[echelonize/GPLU] testing for early abort...\n");
			if (spasm_echelonize_test_completion(A, p, n, U, Uqinv))
				break;
			early_abort_done = 1;
		}
		rows_since_last_pivot += 1;

		/* ensure enough room in L / U for an extra row */
		if (unz + m > U->nzmax)
			spasm_csr_realloc(U, 2 * U->nzmax + m);
		int *Uj = U->j;
		spasm_ZZp *Ux = U->x;
		if (L && lnz + m > L->nzmax)
			spasm_triplet_realloc(L, 2 * L->nzmax + m);
		int *Li = (L != NULL) ? L->i : NULL;
		int *Lj = (L != NULL) ? L->j : NULL;
		spasm_ZZp *Lx = (L != NULL) ? L->x : NULL;

		/* Triangular solve: x * U = A[i] */
		int inew = p[i];
		int i_orig = (p_in != NULL) ? p_in[inew] : inew;
		int top = spasm_sparse_triangular_solve(U, A, inew, xj, x, Uqinv);

		/* Find pivot column; current poor strategy= choose leftmost */
		int jpiv = m ;                 /* column index of best pivot so far. */
		for (int px = top; px < m; px++) {
			int j = xj[px];        /* x[j] is (generically) nonzero */
			if (x[j] == 0)
				continue;
			if (Uqinv[j] < 0) {
				/* non-zero coeff on non-pivotal column --> candidate */
				if (j < jpiv)
					jpiv = j;
			} else if (L != NULL) {
				/* everything under pivotal columns goes into L */
				Li[lnz] = i_orig;
				Lj[lnz] = Uqinv[j];
				Lx[lnz] = x[j];
				lnz += 1;
			}	
		}
		
		if (jpiv == m)
			continue;        /* no pivot found */

		/* add entry entry in L for the pivot */
		if (L != NULL) {
			assert(x[jpiv] != 0);
			Lp[U->n] = i_orig;
			Li[lnz] = i_orig;
			Lj[lnz] = U->n;
			Lx[lnz] = x[jpiv];
			lnz += 1;
		}

		/* store new pivotal row into U */
		Uqinv[jpiv] = U->n;
		// fprintf(stderr, "setting Uqinv[%d] <--- %d\n", jpiv, U->n);

		i64 old_unz = unz;
		Uj[unz] = jpiv;
		Ux[unz] = 1;
		unz += 1;
		// fprintf(stderr, "setting U[%d, %d] <--- 1\n", U->n, jpiv);
		assert(x[jpiv] != 0);
		spasm_ZZp beta = spasm_ZZp_inverse(A->field, x[jpiv]);
		for (int px = top; px < m; px++) {
			int j = xj[px];
			if (x[j] != 0 && Uqinv[j] < 0) {
				Uj[unz] = j;
				Ux[unz] = spasm_ZZp_mul(A->field, beta, x[j]);
				// fprintf(stderr, "setting U[%d, %d] <--- %d\n", U->n, j, Ux[unz]);
				unz += 1;
			}
		}
		U->n += 1;
		Up[U->n] = unz;

		/* reset early abort */
		rows_since_last_pivot = 0;
		early_abort_done = 0;

		if ((i % verbose_step) == 0) {
			fprintf(stderr, "\r[echelonize/GPLU] %d / %d [|U| = %" PRId64 " / |L| = %" PRId64"] -- current density= (%.3f vs %.3f) --- rank >= %d", 
				i, n, unz, lnz, 1.0 * (m - top) / (m), 1.0 * (unz - old_unz) / m, U->n);
			fflush(stderr);
		}
	}
	/* cleanup */
	if (L) {
		L->nz = lnz;
		L->m = U->n;
	}
	fprintf(stderr, "\n");
	free(x);
	free(xj);
}

/*
 * Transfer echelonized rows from (dense) S to (sparse) U
 */
static void update_U_after_rref(int rr, int Sm, const void *S, spasm_datatype datatype, 
	const size_t *Sqinv, const int *q, struct spasm_lu *fact)
{
	struct spasm_csr *U = fact->U;
	int *Uqinv = fact->qinv;
	i64 extra_nnz = ((i64) (1 + Sm - rr)) * rr;     /* maximum size increase */
	i64 unz = spasm_nnz(U);
	fprintf(stderr, "[dense update] enlarging U from %" PRId64 " to %" PRId64 " entries\n", unz, unz + extra_nnz);
	spasm_csr_realloc(U, unz + extra_nnz);
	i64 *Up = U->p;
	int *Uj = U->j;
	spasm_ZZp *Ux = U->x;
	for (i64 i = 0; i < rr; i++) {
		int j = Sqinv[i];   /* column (of S) with the pivot on row i of S; the pivot is implicitly 1 */
		Uj[unz] = q[j];  /* column of A with the pivot */
		Ux[unz] = 1;
		unz += 1;
		Uqinv[q[j]] = U->n;
		for (i64 k = rr; k < Sm; k++) {
			i64 j = Sqinv[k];
			spasm_ZZp x = spasm_datatype_read(S, i * Sm + k, datatype);
			if (x == 0)
				continue;   /* don't store zero */
			Uj[unz] = q[j];
			Ux[unz] = x;  // reduce?
			unz += 1;
		}
		U->n += 1;
		Up[U->n] = unz;
	}
	assert(unz == spasm_nnz(U));
}

/*
 * Transfer dense LU factorization to fact
 */
static void update_fact_after_LU(int n, int Sm, int r, const void *S, spasm_datatype datatype, 
	const size_t *Sp, const size_t *Sqinv, const int *q, const int *p_in, i64 lnz_before, 
	bool complete, bool *pivotal, struct spasm_lu *fact)
{
	struct spasm_csr *U = fact->U;
	struct spasm_triplet *L = fact->Ltmp;
	int *Uqinv = fact->qinv;
	int *Lp = fact->p;
	i64 extra_unz = ((i64) (1 + 2*Sm - r)) * r;     /* maximum size increase */
	i64 extra_lnz = ((i64) (2*n - r + 1)) * r / 2;
	i64 unz = spasm_nnz(U);
	i64 lnz = L->nz;
	spasm_csr_realloc(U, unz + extra_unz);
	spasm_triplet_realloc(L, lnz + extra_lnz);
	i64 *Up = U->p;
	int *Uj = U->j;
	spasm_ZZp *Ux = U->x;
	int *Li = L->i;
	int *Lj = L->j;
	spasm_ZZp *Lx = L->x;
	
	/* build L */
	if (!complete) {
		for (i64 i = 0; i < r; i++) {   /* mark pivotal rows */
			int pi = Sp[i];
			int iorig = (p_in != NULL) ? p_in[pi] : pi;
			pivotal[iorig] = 1;
		}

		/* stack L (ignore non-pivotal rows) */
		lnz = lnz_before;
		for (i64 px = lnz_before; px < L->nz; px++) {
			int i = Li[px];
			if (!pivotal[i])
				continue;
			int j = Lj[px];
			spasm_ZZp x = Lx[px];
			Li[lnz] = i;
			Lj[lnz] = j;
			Lx[lnz] = x;
			lnz += 1;
		}
		fprintf(stderr, "L : %" PRId64 " --> %" PRId64 " --> %" PRId64 " ---> ", lnz_before, L->nz, lnz);
	}

	/* add new entries from S */
	for (i64 i = 0; i < (complete ? n : r); i++) {
		int pi = Sp[i];
		int iorig = (p_in != NULL) ? p_in[pi] : pi;
		for (i64 j = 0; j < spasm_min(i + 1, r); j++) {
			spasm_ZZp Mij = spasm_datatype_read(S, i  * Sm + j, datatype);
			if (Mij == 0)
				continue;
			Li[lnz] = iorig;
			Lj[lnz] = U->n + j;
			Lx[lnz] = Mij;
			lnz += 1;
		}
		if (i < r)   /* register pivot */
			Lp[U->n + i] = iorig;
	}
	L->nz = lnz;
	fprintf(stderr, "%" PRId64 "\n", lnz);

	/* fill U */
	for (i64 i = 0; i < r; i++) {
		/* implicit 1 in U */
		int j = Sqinv[i];
		int jorig = q[j];
		Uj[unz] = jorig;
		Ux[unz] = 1;
		unz += 1;
		/* register pivot */
		Uqinv[jorig] = U->n;
		for (i64 j = i+1; j < Sm; j++) {
			int jnew = Sqinv[j];
			int jorig = q[jnew];
			spasm_ZZp x = spasm_datatype_read(S, i * Sm + j, datatype);
			Uj[unz] = jorig;
			Ux[unz] = x;
			unz += 1;
		}
		U->n += 1;
		Up[U->n] = unz;
	}
}

static void echelonize_dense_lowrank(const struct spasm_csr *A, const int *p, int n, struct spasm_lu *fact, struct echelonize_opts *opts)
{
	assert(opts->dense_block_size > 0);
	struct spasm_csr *U = fact->U;
	int *Uqinv = fact->qinv;
	int m = A->m;
	int Sm = m - U->n;
	i64 prime = spasm_get_prime(A);
	spasm_datatype datatype = spasm_datatype_choose(prime);

	i64 size_S = (i64) opts->dense_block_size * (i64) Sm * spasm_datatype_size(datatype);
	void *S = spasm_malloc(size_S);
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sp = spasm_malloc(Sm * sizeof(*Sp));       /* for FFPACK */
	double start = spasm_wtime();
	int old_un = U->n;
	int round = 0;
	fprintf(stderr, "[echelonize/dense/low-rank] processing dense schur complement of dimension %d x %d; block size=%d, type %s\n", 
		n, Sm, opts->dense_block_size, spasm_datatype_name(datatype));
	
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
		int Sn = spasm_min(rank_ub, opts->dense_block_size);
		if (Sn <= 0)
			break;		
		fprintf(stderr, "[echelonize/dense/low-rank] Round %d. Weight %d. Processing chunk (%d x %d), |U| = %"PRId64"\n", 
			round, w, Sn, Sm, spasm_nnz(U));
		spasm_schur_dense_randomized(A, p, n, U, Uqinv, S, datatype, q, Sn, w);
		int rr = spasm_ffpack_rref(prime, Sn, Sm, S, Sm, datatype, Sp);

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
		update_U_after_rref(rr, Sm, S, datatype, Sp, q, fact);
		n -= rr;
		Sm -= rr;
		rank_ub -= rr;
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
static void echelonize_dense(const struct spasm_csr *A, const int *p, int n, const int *p_in, struct spasm_lu *fact, struct echelonize_opts *opts)
{
	assert(opts->dense_block_size > 0);
	struct spasm_csr *U = fact->U;
	int m = A->m;
	int Sm = m - U->n;
	i64 prime = spasm_get_prime(A);
	spasm_datatype datatype = spasm_datatype_choose(prime);

	void *S = spasm_malloc((i64) opts->dense_block_size * Sm * spasm_datatype_size(datatype));
	int *p_out = spasm_malloc(opts->dense_block_size * sizeof(*p_out));
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sqinv = spasm_malloc(Sm * sizeof(*Sqinv));                   /* for FFPACK */
	size_t *Sp = spasm_malloc(opts->dense_block_size * sizeof(*Sp));     /* for FFPACK / LU only */
	bool *pivotal = spasm_malloc(A->n * sizeof(*pivotal));
	for (int i = 0; i < A->n; i++)
		pivotal[i] = 0;

	int processed = 0;
	double start = spasm_wtime();
	int old_un = U->n;
	int round = 0;
	fprintf(stderr, "[echelonize/dense] processing dense schur complement of dimension %d x %d; block size=%d, type %s\n", 
		n, Sm, opts->dense_block_size, spasm_datatype_name(datatype));
	bool lowrank_mode = 0;
	int rank_ub = spasm_min(A->n - U->n, A->m - U->n);

	for (;;) {
		/* compute a chunk of the schur complement, then echelonize with FFPACK */
		int Sn = spasm_min(opts->dense_block_size, n - processed);
		if (Sn <= 0)
			break;
		
		fprintf(stderr, "[echelonize/dense] Round %d. processing S[%d:%d] (%d x %d)\n", round, processed, processed + Sn, Sn, Sm);	

		i64 lnz_before = (opts->L) ? fact->Ltmp->nz : -1;
		spasm_schur_dense(A, p, Sn, p_in, fact, S, datatype, q, p_out);

		int rr;
		if (opts->L) {
			rr = spasm_ffpack_LU(prime, Sn, Sm, S, Sm, datatype, Sp, Sqinv);
			update_fact_after_LU(Sn, Sm, rr, S, datatype, Sp, Sqinv, q, p_out, lnz_before, opts->complete, pivotal, fact);
		} else {
			rr = spasm_ffpack_rref(prime, Sn, Sm, S, Sm, datatype, Sqinv);
			update_U_after_rref(rr, Sm, S, datatype, Sqinv, q, fact);
		}

		// TODO: test completion and allow early abort

		/* move on to the next chunk */
		round += 1;
		processed += Sn;
		p += Sn;
		Sm = m - U->n;
		rank_ub = spasm_min(A->n - U->n, A->m - U->n);
		fprintf(stderr, "[echelonize/dense] found %d new pivots\n", rr);

		/* 
		 * switch to low-rank mode if yield drops too much 
		 * This will early abort if a full factorization is not needed
		 */
		if (opts->enable_tall_and_skinny && (rr < opts->low_rank_ratio * Sn)) {
			lowrank_mode = 1;
			break;
		}
	}
	free(S);
	free(q);
	free(Sqinv);
	free(Sp);
	free(p_out);
	free(pivotal);
	if (rank_ub > 0 && n - processed > 0 && lowrank_mode) {
		fprintf(stderr, "[echelonize/dense] Too few pivots; switching to low-rank mode\n");
		echelonize_dense_lowrank(A, p, n - processed, fact, opts);
	} else {
		fprintf(stderr, "[echelonize/dense] completed in %.1fs. %d new pivots found\n", spasm_wtime() - start, U->n - old_un);
	}
}


/*
 * (main entry point)
 * Returns the row echelon form of A. 
 * Initializes Uqinv (must be preallocated of size m [==#columns of A]).
 * Modifies A (permutes entries in rows)
 * FIXME potential memleak (>= 1 rounds then status == 1...)
 */
struct spasm_lu * spasm_echelonize(const struct spasm_csr *A, struct echelonize_opts *opts)
{
	struct echelonize_opts default_opts;
	if (opts == NULL) {
		fprintf(stderr, "[echelonize] using default settings\n");
		opts = &default_opts;
		spasm_echelonize_init_opts(opts);
	}
	int n = A->n;
	int m = A->m;
	i64 prime = spasm_get_prime(A);
	fprintf(stderr, "[echelonize] Start on %d x %d matrix with %" PRId64 " nnz\n", n, m, spasm_nnz(A));
	
	/* options sanity check */
	if (opts->complete)
		opts->L = 1;
	if (opts->L)
		opts->enable_tall_and_skinny = 0;   // for now

	/* allocate result */
	struct spasm_csr *U = spasm_csr_alloc(n, m, spasm_nnz(A), prime, true);
	int *Uqinv = spasm_malloc(m * sizeof(*Uqinv));
	U->n = 0;
	for (int j = 0; j < m; j++)
		Uqinv[j] = -1;
	
	struct spasm_triplet *L = NULL;
	int *Lp = NULL;
	if (opts->L) {
		L = spasm_triplet_alloc(n, n, spasm_nnz(A), prime, true);
		Lp = spasm_malloc(n * sizeof(*Lp));
		for (int j = 0; j < n; j++)
			Lp[j] = -1;
		assert(L->x != NULL);
	}
	
	struct spasm_lu *fact = spasm_malloc(sizeof(*fact));
	fact->L = NULL;
	fact->p = Lp;
	fact->U = U;
	fact->qinv = Uqinv;
	fact->Ltmp = L;

	/* local stuff */
	int *p = spasm_malloc(n * sizeof(*p)); /* pivotal rows come first in P*A */
	double start = spasm_wtime();
	double density = (double) spasm_nnz(A) / n / m;
	int npiv = 0;
	int status = 0;  /* 0 == max_round reached; 1 == full rank reached; 2 == early abort */
	int *p_in = NULL;

	int round;
	for (round = 0; round < opts->max_round; round++) {
		/* decide whether to move on to the next iteration */
		if (spasm_nnz(A) == 0) {
			fprintf(stderr, "[echelonize] empty matrix\n");
			status = 1;
			break;
		}
		
		fprintf(stderr, "[echelonize] round %d\n", round);
		npiv = spasm_pivots_extract_structural(A, p_in, fact, p, opts);

		if (npiv < opts->min_pivot_proportion * spasm_min(n, m - U->n)) {
			fprintf(stderr, "[echelonize] not enough pivots found; stopping\n");
			status = 2;
			break;     /* not enough pivots found */
		}		
		// if (density > opts->sparsity_threshold && aspect_ratio > opts->tall_and_skinny_ratio) {
		// 	fprintf(stderr, "Schur complement is dense, tall and skinny (#rows / #cols = %.1f)\n", aspect_ratio);
		// 	break;
		// }
		density = spasm_schur_estimate_density(A, p + npiv, n - npiv, U, Uqinv, 100);
		if (density > opts->sparsity_threshold) {
			fprintf(stderr, "[echelonize] Schur complement is dense (estimated %.2f%%)\n", 100 * density);
			status = 2;
			break;
		}

		/* compute the next schur complement */
		i64 nnz = (density * (n - npiv)) * (m - U->n);
		char tmp[8];
		spasm_human_format(sizeof(int) * (n - npiv + nnz) + sizeof(spasm_ZZp) * nnz, tmp);
		fprintf(stderr, "Schur complement is %d x %d, estimated density : %.2f (%s byte)\n", n - npiv, m - U->n, density, tmp);
		int *p_out = spasm_malloc((n - npiv) * sizeof(*p_out));
		struct spasm_csr *S = spasm_schur(A, p + npiv, n - npiv, fact, density, L, p_in, p_out);
		if (round > 0)
			spasm_csr_free((struct spasm_csr *) A);       /* discard const, only if it is not the input argument */
		A = S;
		n = n - npiv;
		free(p_in);
		p_in = p_out;
	}
	/*
	 * status == 0. Exit because opts->max_round reached. Just factor A.
	 * status == 1. Exit because A == 0. Nothing more to do.
	 * status == 2. Some pivots found (U/L updated), but schur complement not computed (too few pivots / too dense).
	 */

	if (status == 0) {
		npiv = 0;
		for (int i = 0; i < n; i++)
			p[i] = i;
	}
	if (status == 1)
		goto cleanup;  /* nothing else to do */

	/* finish */
	if (!opts->enable_tall_and_skinny)
		fprintf(stderr, "[echelonize] dense low-rank mode disabled\n");
	if (!opts->enable_dense)
		fprintf(stderr, "[echelonize] regular dense mode disabled\n");
	if (!opts->enable_GPLU)
		fprintf(stderr, "[echelonize] GPLU mode disabled\n");
	
	double aspect_ratio = (double) (n - npiv) / (m - U->n);
	fprintf(stderr, "[echelonize] finishing; density = %.3f; aspect ratio = %.1f\n", density, aspect_ratio);
	if (opts->enable_tall_and_skinny && aspect_ratio > opts->tall_and_skinny_ratio)
		echelonize_dense_lowrank(A, p + npiv, n - npiv, fact, opts);
	else if (opts->enable_dense && density > opts->sparsity_threshold)
		echelonize_dense(A, p + npiv, n - npiv, p_in, fact, opts);
	else if (opts->enable_GPLU)
		echelonize_GPLU(A, p + npiv, n - npiv, p_in, fact, opts);
	else
		fprintf(stderr, "[echelonize] Cannot finish (no valid method enabled). Incomplete echelonization returned\n");

cleanup:
	free(p);
	free(p_in);
	fprintf(stderr, "[echelonize] Done in %.1fs. Rank %d, %" PRId64 " nz in basis\n", spasm_wtime() - start, U->n, spasm_nnz(U));
	spasm_csr_resize(U, U->n, m);
	spasm_csr_realloc(U, -1);
	if (round > 0)
		spasm_csr_free((struct spasm_csr *) A);
	if (opts->L) {
		L->m = U->n; 
		fact->p = spasm_realloc(Lp, U->n * sizeof(*Lp));
		fact->L = spasm_compress(L);
		spasm_triplet_free(L);
		fact->Ltmp = NULL;
		fact->complete = opts->complete;
	}
	fact->r = U->n;
	return fact;
}