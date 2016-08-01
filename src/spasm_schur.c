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
 * p[0] ... p[n_pivots-1] of input matrix A. The pivots must be the entries
 * on the lines. This returns a sparse representation of S. The pivots must
 * be unitary.
 */
spasm *spasm_schur(spasm * A, const int *p, const int *qinv, const int npiv) {
	spasm *S;
	int k, *Sp, *Sj, Sn, Sm, m, n, snz, *xj, top, *q, verbose_step;
	spasm_GFp *Sx, *x;

	/* check inputs */
	assert(A != NULL);
	n = A->n;
	m = A->m;
	assert(n >= npiv);
	assert(m >= npiv);

	/* Get Workspace */
	Sn = n - npiv;
	Sm = m - npiv;
	snz = 4 * (Sn + Sm);	/* educated guess */
	S = spasm_csr_alloc(Sn, Sm, snz, A->prime, SPASM_WITH_NUMERICAL_VALUES);

	x = spasm_malloc(m * sizeof(spasm_GFp));
	xj = spasm_malloc(3 * m * sizeof(int));
	spasm_vector_zero(xj, 3 * m);

	verbose_step = spasm_max(1, n / 1000);
	Sp = S->p;
	Sj = S->j;
	Sx = S->x;

	/*
	 * q sends the non-pivotal columns of A to the columns of S. It is
	 * not the inverse of qinv...
	 */
	q = spasm_malloc(m * sizeof(int));
	k = 0;
	for (int j = 0; j < m; j++)
		q[j] = (qinv[j] < 0) ? k++ : -1;

	snz = 0;		/* non-zero in S */
	Sn = 0;			/* rows in S */

	fprintf(stderr, "Starting Schur complement computation...\n");
	for (int i = npiv; i < n; i++) {
		int inew = p[i];

		/* triangular solve */
		top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);

		/* not enough room in S ? realloc twice the size */
		if (snz + Sm > S->nzmax) {
			spasm_csr_realloc(S, 2 * S->nzmax + Sm);
			Sj = S->j;
			Sx = S->x;
		}
		Sp[Sn] = snz;	/* S[i] starts here */
		for (int px = top; px < m; px++) {
			int j = xj[px];

			if (x[j] == 0)
				continue;

			/* save non-zero, non-pivot coefficients in S */
			if (q[j] >= 0) {
				Sj[snz] = q[j];
				Sx[snz] = x[j];
				snz++;
			}
		}
		Sn++;

		if ((i % verbose_step) == 0) {
			fprintf(stderr, "\rSchur : %d / %d [S=%d * %d, %d NNZ] -- current density= (%.3f)", i, n, Sn, Sm, snz, 1.0 * snz / (1.0 * Sm * Sn));
			fflush(stderr);
		}
	}
	/* finalize S */
	fprintf(stderr, "\n");
	Sp[S->n] = snz;
	spasm_csr_realloc(S, -1);

	/* free extra workspace */
	free(q);
	free(x);
	free(xj);

	return S;
}


/** Samples R rows at random in the schur complement of A w.r.t. the pivots in p[0:n_pivots],
* and return the number that are non-zero (these rows of A are linearly independent from the pivots).
* The pivots must be unitary.
*/
double spasm_schur_probe_density(spasm * A, const int *p, const int *qinv, const int npiv, const int R) {
	int m, n, *xj, top, nnz;
	spasm_GFp *x;

	/* check inputs */
	m = A->m;
	n = A->n;

	/* Get Workspace */
	x = spasm_malloc(m * sizeof(spasm_GFp));
	xj = spasm_malloc(3 * m * sizeof(int));
	spasm_vector_zero(xj, 3 * m);

	nnz = 0;
	for (int i = 0; i < R; i++) {
		/* pick a random row in S, check if non-zero */
		int inew = p[npiv + (rand() % (n - npiv))];
		top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);
		for (int px = top; px < m; px++) {
			int j = xj[px];
			if (qinv[j] < 0 && x[j] != 0)
				nnz++;
		}
	}

	/* free extra workspace */
	free(x);
	free(xj);

	return ((double)nnz) / (m - npiv) / R;
}

/*
 * computes the rank of the schur complement, but not the schur complement
 * itself. The pivots must be unitary.
 */
int spasm_schur_rank(spasm * A, const int *p, const int *qinv, const int npiv) {
	int Sn, Sm, m, n, k, r, prime, new, step, nbad, ngood, abort;
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
	Sn = n - npiv;
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
	nbad = 0;
	ngood = 0;
	abort = 0;

#pragma omp parallel
	{
		spasm_GFp * x = spasm_malloc(m * sizeof(spasm_GFp));
		spasm_GFp * y = spasm_malloc(Sm * sizeof(spasm_GFp));

		while (!abort) {
			/* l'essentiel du temps se passe avec step == 1 */

			/* make sure I have the latest value of step */
			#pragma omp flush(step)

			/* random linear combination */
			spasm_vector_zero(x, m);
			if (step >= n)	/* on prend tout ! */
				for (int i = npiv; i < n; i++) {
					int inew = p[i];
					spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, x, prime);
				}
			else		/* on prend un sous-ensemble aléatoire */
				for (int i = 0; i < step; i++) {
					int inew = p[npiv + (rand() % (n - npiv))];
					spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], 1 + (rand() % (prime - 1)), x, prime);
				}

			spasm_eliminate_sparse_pivots(A, npiv, p, x);

			/* Gather the solution in y */
			for (int j = 0; j < Sm; j++) {
				y[j] = x[q[j]];
			}

			new = spasm_dense_LU_process(U, y);

#pragma omp critical(schur_rank)
			if (new) {
				r++;
				ngood++;
				if (ngood == 16) {
					if (step > 1) {
						nbad = 0;
					}
					step = spasm_max(step / 2, 1);
					ngood = 0;
				}
			} else {
				nbad++;
				if (nbad == 3) {
					if (step > n) {
						abort = 1;	/* c'est fini */
					}
					step *= 2;
					nbad = 0;
					ngood = 0;
				}
			}

			#pragma omp atomic update
			k++;
			fprintf(stderr, "\rSchur : %d [step = %d] -- current rank = %d / bad = %d, good = %d   ", k, step, r, nbad, ngood);
			fflush(stderr);

			#pragma omp flush(abort)
		}
		free(x);
		free(y);
	}
	fprintf(stderr, "\n[schur/rank] Time: %.1fs\n", spasm_wtime() - start);
	free(q);
	spasm_dense_LU_free(U);
	return r;
}
