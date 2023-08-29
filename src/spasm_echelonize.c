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
static const double sparsity_threshold = 0.05;

/* for debugging purposes */
void validate_pivots(const spasm *A, int npiv, const int *qinv, const int *p)
{
	int *Ap = A->p;
	int *Aj = A->j;
	spasm_GFp *Ax = A->x;
	for (int k = 0; k < npiv; k++) {
		int i = (p != NULL) ? p[k] : k;
		int j = Aj[Ap[i]];
		assert(Ax[Ap[i]] == 1);
		assert(qinv[j] == i); 
	}
}



spasm* spasm_echelonize(spasm *A, int *Uqinv)
{
	int n = A->n;
	int m = A->m;
	int prime = A->prime;
	
	int *qinv = spasm_malloc(n * sizeof(int));  /* for pivot search */
	int *p = spasm_malloc(n * sizeof(int));

	spasm *U = spasm_csr_alloc(n, m, spasm_nnz(A), prime, SPASM_WITH_NUMERICAL_VALUES);
	U->n = 0;
	int Unz = 0;      /* #entries in U */
	for (int j = 0; j < m; j++)
		Uqinv[j] = -1;

	double start = spasm_wtime();
	int round = 0;
	double density = -1;
	for (;;) {
	 	int *Ap = A->p;
		int *Aj = A->j;
		spasm_GFp *Ax = A->x;

		/* invariant: A does not have coefficients on pivotal columns */
		for (int i = 0; i < n; i++)
			for (int px = Ap[i]; px < Ap[i + 1]; px++) {
				int j = Aj[px];
				assert(Uqinv[j] < 0);
			}

		/* find structural pivots in A */
		int npiv = spasm_find_pivots(A, p, qinv);

		if (npiv < 0.01 * spasm_min(n, m)) {
			fprintf(stderr, "not enough pivots found\n");
			/* TODO: echelonize A, take all the new pivots into U */
			break;     /* not enough pivots found */
		}

		/* compute total pivot nnz and reallocate U if necessary */
		int pivot_nnz = 0;
		for (int k = 0; k < npiv; k++) {
			int i = p[k];
			pivot_nnz += spasm_row_weight(A, i);
		}
		if (spasm_nnz(U) + pivot_nnz >= U->nzmax)
			spasm_csr_realloc(U, spasm_nnz(U) + pivot_nnz);

		/* copy pivotal rows to U and make them unitary; update Uqinv */
		int *Up = U->p;
		int *Uj = U->j;
		spasm_GFp *Ux = U->x;

		for (int k = 0; k < npiv; k++) {
			int i = p[k];
			int j = Aj[Ap[i]];
			assert(Uqinv[j] < 0);
			assert(qinv[j] == i);
			Uqinv[j] = U->n;
			for (int px = Ap[i]; px < Ap[i + 1]; px++) {
				Uj[Unz] = Aj[px];
				Ux[Unz] = Ax[px];
				Unz += 1;
			}
			U->n += 1;
			Up[U->n] = Unz;
		}
		spasm_make_pivots_unitary(U, SPASM_IDENTITY_PERMUTATION, U->n); /* reprocess previous pivots; shame */
		validate_pivots(U, U->n, Uqinv, NULL);

		density = estimate_schur_density(A, npiv, p, U, Uqinv, 100);
		if (density > sparsity_threshold) {
			fprintf(stderr, "Schur complement is dense (esitmated %.2f%%)\n", 100 * density);
			/* TODO: compute the dense schur complement and echelonize it */
			break;     /* schur complement is not sparse enough */
		}

		int64_t nnz = (density * (n - npiv)) * (m - npiv);
		char tmp[6];
		spasm_human_format(sizeof(int) * (n - npiv + nnz) + sizeof(spasm_GFp) * nnz, tmp);
		fprintf(stderr, "round %d. Schur complement is %d x %d, estimated density : %.2f (%s byte)\n", round, n - npiv, m - npiv, density, tmp);

		/* compute the schur complement */
		spasm *S = spasm_schur(A, p, npiv, U, Uqinv, density, SPASM_DISCARD_L, NULL);
		A = S;
		n = A->n;
		round += 1;
	}

	/* Final step: echelonize A, rows p[npiv:n] */

	/* sparse schur complement : GPLU */
	if (1 || density < sparsity_threshold) {
		spasm_lu *LU = spasm_LU(A, p, SPASM_DISCARD_L);
		spasm *UU = LU->U; 
		int *UUp = UU->p;
		int *UUj = UU->j;
		spasm_GFp *UUx = UU->x;
		int *qinv = LU->qinv;
		// concatenate UU to U;
		spasm_csr_realloc(U, spasm_nnz(U) + spasm_nnz(UU));
		int *Up = U->p;
		int *Uj = U->j;
		spasm_GFp *Ux = U->x;
		for (int i = 0; i < UU->n; i++) {
			int j = UUj[UUp[i]];
			assert(qinv[j] == i);
			Uqinv[j] = U->n;
			for (int px = UUp[i]; px < UUp[i + 1]; px++) {
				Uj[Unz] = UUj[px];
				Ux[Unz] = UUx[px];
				Unz += 1;
			}
			U->n += 1;
			Up[U->n] = Unz;   /* finalization */
		}
		spasm_free_LU(LU);
	}// else {
//		/* dense schur complement */
//		int r = spasm_schur_rank(A, p, qinv, npiv);
//		fprintf(stderr, "rank = %d + %d\n", npiv, r);
	//}
	return U;
}