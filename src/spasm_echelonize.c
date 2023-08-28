#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

#include "spasm.h"

/** NOT DRY WRT spasm_schur
 * Samples R rows at random in the schur complement of A w.r.t. U, and return the average density.
 * rows p[0:npiv] of A are pivotal and must be excluded.
 * qinv locates the pivots in U.
 */
double estimate_schur_density(const spasm * A, int npiv, const int *p, const spasm *U, const int *qinv, int R)
{
	int m = A->m;
	int n = A->n;
	int nnz = 0;
	assert(npiv < m);

	#pragma omp parallel
	{
		/* per-thread scratch space */
		spasm_GFp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * m * sizeof(*xj));
		spasm_vector_zero(xj, 3 * m);

		#pragma omp for reduction(+:nnz) schedule(dynamic)
		for (int i = 0; i < R; i++) {
			/* pick a random non-pivotal row in A */
			int inew = p[npiv + (rand() % (n - npiv))];
			int top = spasm_sparse_forward_solve(U, A, inew, xj, x, qinv);
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (qinv[j] < 0 && x[j] != 0)
					nnz += 1;
			}
		}
		free(x);
		free(xj);
	}
	return ((double) nnz) / (m - npiv) / R;
}

/*
 * Returns the row echelon form of A. Initializes qinv (must be of size m, #columns(A)).
 * Modifies A (permutes entries only)
 */
spasm* spasm_echelonize(spasm *A, int *Uqinv)
{
	int n = A->n;
	int m = A->m;
	int *Ap = A->p;
	int *Aj = A->j;
	spasm_GFp *Ax = A->x;
	int prime = A->prime;
	
	int *qinv = spasm_malloc(n * sizeof(int));  /* for pivot search */
	int *p = spasm_malloc(n * sizeof(int));

	spasm *U = spasm_csr_alloc(n, m, spasm_nnz(A), prime, SPASM_WITH_NUMERICAL_VALUES);
	int Unz = 0;      /* #entries in U */
	int Un = 0;       /* #rows in U */
	for (int j = 0; j < m; j++)
		Uqinv[j] = -1;

	double start = spasm_wtime();	
	int round = 0;
	for (;;) {
		/* find structural pivots in A */
		int npiv = spasm_find_pivots(A, p, qinv);

		/* compute total pivot nnz and reallocate U if necessary */
		int pivot_nnz = 0;
		for (int k = 0; k < npiv; k++) {
			int i = p[k];
			pivot_nnz += spasm_row_weigh(A[i]);
		}
		if (U->nz + pivot_nnz >= U->nzmax)
			spasm_csr_realloc(U, U->nz + pivot_nnz);

		/* copy pivotal rows to U and make them unitary; update Uqinv */
		int *Up = U->p;
		int *Uj = U->j;
		spasm_GFp *Ux = U->x;
		for (int k = 0; k < npiv; k++) {
			Up[Un] = Unz;
			int i = p[k];
			assert(qinv[Ap[i]] == i);
			qinv[Ap[i]] = k;
			for (int px = Ap[i]; px < Ap[i + 1]; px++) {
				Uj[Unz] = Aj[px];
				Ux[Unz] = Ax[px];
				Unz += 1;
			}
			Un += 1;
		}
		Up[Un] = nz;   /* finalization */
		spasm_make_pivots_unitary(U, SPASM_IDENTITY_PERMUTATION, Un); /* reprocess previous pivots; shame */
	
		/* abort if the schur complement is not sparse enough */
		double density = estimate_schur_density(A, npiv, p, U, Uqinv, 100);
		if (density > sparsity_threshold)
			break;

		if (npiv < 0.01 * spasm_min(n, m))
			break;     /* abort if not enough pivots are found */

		int64_t nnz = (density * (n - npiv)) * (m - npiv);
		char tmp[6];
		spasm_human_format(sizeof(int) * (n - npiv + nnz) + sizeof(spasm_GFp) * nnz, tmp);
		fprintf(stderr, "round %d. Schur complement is %d x %d, estimated density : %.2f (%s byte)\n", round, n - npiv, m - npiv, density, tmp);

		/* compute the schur complement */
		spasm *S = spasm_schur(A, p, npiv, density, 0, NULL); ///////// FIXME must read U
		A = S;
		n = A->n;
	}

	/* Final step: echelonize A, except rows p[0:npiv] */

	/* sparse schur complement : GPLU */
	if (gplu_final || (!dense_final && density < sparsity_threshold)) {
		spasm_lu *LU = spasm_LU(A, p, SPASM_DISCARD_L);
		rank += LU->U->n;
		spasm_free_LU(LU);
	} else {
		/* dense schur complement */
		int r = spasm_schur_rank(A, p, qinv, npiv);
		fprintf(stderr, "rank = %d + %d\n", npiv, r);
		rank += npiv + r;
	}
