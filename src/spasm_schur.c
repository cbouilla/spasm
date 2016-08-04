/* indent -nfbs -nip -npsl -di0 spasm_schur.c  */
#include <assert.h>
#include "spasm.h"

/* make pivotal rows of A unitary */
void spasm_make_pivots_unitary(spasm * A, const int *p, const int npiv) {
	int prime = A->prime;
	int *Ap = A->p;
	spasm_GFp *Ax = A->x;

#pragma omp parallel for
	for (int i = 0; i < npiv; i++) {
		int inew;
		spasm_GFp diag, alpha;

		inew = p[i];
		diag = Ax[Ap[inew]];
		if (diag == 1)
			continue;

		alpha = spasm_GFp_inverse(diag, prime);
		for (int px = Ap[inew]; px < Ap[inew + 1]; px++)
			Ax[px] = (alpha * Ax[px]) % prime;
	}
}

/*
 * Computes the Schur complement, by eliminating the pivots located on rows
 * p[0] ... p[n_pivots-1] of input matrix A. The pivots must be the first entries
 * on the lines. This returns a sparse representation of S. The pivots must
 * be unitary.
 */
spasm *spasm_schur(spasm * A, const int *p, const int *qinv, const int npiv) {
	int k, snz, *q, writing;
	double start;
	spasm *S;

	const int n = A->n;
	const int m = A->m;
	const int Sm = m - npiv;
	const int Sn = n - npiv;
	const int verbose_step = spasm_max(1, n / 1000);
	
	/*
	 * q sends the non-pivotal columns of A to the columns of S. It is
	 * not the inverse of qinv...
	 */
	q = spasm_malloc(m * sizeof(int));
	k = 0;
	for (int j = 0; j < m; j++)
		q[j] = (qinv[j] < 0) ? k++ : -1;

	S = spasm_csr_alloc(Sn, Sm, Sn + Sm, A->prime, SPASM_WITH_NUMERICAL_VALUES);
	int *Sp = S->p;
	int *Sj = S->j;
	spasm_GFp *Sx = S->x;
	snz = 0;
	k = 0;
	writing = 0;
	start = spasm_wtime();

#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
		int *xj = spasm_malloc(3 * m * sizeof(int));
		spasm_vector_zero(xj, 3 * m);		
		int tid = spasm_get_thread_num();
		int row_snz, row_k, row_px;

		#pragma omp for schedule(dynamic, verbose_step)
		for (int i = npiv; i < n; i++) {
			const int inew = p[i];
			const int top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);

			row_snz = 0;
			for (int px = top; px < m; px++) {
				const int j = xj[px];
				if ((q[j] >= 0) && (x[j] != 0))
					row_snz++;
			}

			#pragma omp critical(schur_complement)
			{
				/* not enough room in S ? realloc twice the size */
				if (snz + Sm > S->nzmax) {
					/* TODO wait until other threads stop writing into it */
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
			
			/* write the new row in S */
			Sp[row_k] = row_px;
			for (int px = top; px < m; px++) {
				const int j = xj[px];
				if ((q[j] >= 0) && (x[j] != 0)) {
					Sj[row_px] = q[j];
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
	free(q);
	return S;
}


/** Samples R rows at random in the schur complement of A w.r.t. the pivots in p[0:n_pivots],
* and return the number that are non-zero (these rows of A are linearly independent from the pivots).
* The pivots must be unitary.
*/
double spasm_schur_probe_density(spasm * A, const int *p, const int *qinv, const int npiv, const int R) {
	int nnz;
	
	const int m = A->m;
	const int n = A->n;
	nnz = 0;

#pragma omp parallel 
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
		int *xj = spasm_malloc(3 * m * sizeof(int));
		spasm_vector_zero(xj, 3 * m);

#pragma omp for reduction(+:nnz) schedule(dynamic)
		for (int i = 0; i < R; i++) {
			/* pick a random row in S, check if non-zero */
			int inew = p[npiv + (rand() % (n - npiv))];
			int top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);
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
int spasm_schur_rank(spasm * A, const int *p, const int *qinv, const int npiv) {
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
