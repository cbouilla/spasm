#include <assert.h>

#include "spasm.h"

/** produces a matrix U such that:
    *) rowspan(U) == rowspan(A) 
    *) there exists a column permutation Q such that U.Q is upper-trapezoidal 

    THIS DESTROYS A.
*/
spasm * spasm_echelonize(spasm *A, int n_iteration)
{
	if (n_iteration < 0)
		n_iteration = 3;
	
	int n = A->n;
	int m = A->m;

	/* allocate room for U */
	spasm *U = spasm_csr_alloc(n, m, spasm_nnz(A), A->prime, SPASM_WITH_NUMERICAL_VALUES);
	
	int *p = spasm_malloc(n * sizeof(int));
	int *qinv = spasm_malloc(m * sizeof(int));
	int *Up = U->p;
	int *Uj = U->j;
	spasm_GFp *Ux = U->x;
	int unz = 0;
	int u_n = 0;

	double start_time = spasm_wtime();

	for (int iteration = 0; iteration < n_iteration; iteration++) {
		/* TODO : do something if A is empty... */
		
		/* find new pivots */
		int npiv = spasm_find_pivots(A, p, qinv);
		spasm_make_pivots_unitary(A, p, npiv);
	
		/* copy them into U */
		int *Ap = A->p;
		int *Aj = A->j;
		int *Ax = A->x;	
		for (int i = 0; i < npiv; i++) {
			Up[u_n] = unz;
			int I = p[i];
			/* not enough room in U ? realloc twice the size */
			if (unz + m > U->nzmax) {
				spasm_csr_realloc(U, 2 * U->nzmax + m);
				Uj = U->j;
				Ux = U->x;
			}
			for (int px = Ap[I]; px < Ap[I + 1]; px++) {
				Uj[unz] = Aj[px];
				Ux[unz] = Ax[px];
				unz++;
			}
			u_n++;
		}
		Up[u_n] = unz;

		/* compute schur complement, update matrix */
		spasm *B = spasm_schur(A, p, npiv, -1, 0, NULL);
		spasm_csr_free(A);
		A = B;
	}

	/* final step : GPLU on the remainder */
	int npiv = spasm_find_pivots(A, p, qinv);
	spasm_make_pivots_unitary(A, p, npiv);
	spasm_lu *LU = spasm_GPLU(A, p, SPASM_DISCARD_L);
	
	/* append last pivots to U */
	spasm_csr_free(A);
	A = LU->U;
	
	int *Ap = A->p;
	int *Aj = A->j;
	int *Ax = A->x;	
	for (int i = 0; i < A->n; i++) {
		/* not enough room in U ? realloc twice the size */
		if (unz + m > U->nzmax) {
			spasm_csr_realloc(U, 2 * U->nzmax + m);
			Uj = U->j;
			Ux = U->x;
		}

		Up[u_n] = unz;
		for (int px = Ap[i]; px < Ap[i + 1]; px++) {
			Uj[unz] = Aj[px];
			Ux[unz] = Ax[px];
			unz++;
		}
		u_n++;
	}
	Up[u_n] = unz;
	U->n = u_n;
	spasm_free_LU(LU);

	char nnz[6];
	spasm_human_format(spasm_nnz(U), nnz);
	fprintf(stderr, "[echelonize] done in %.3f s. NNZ(U) = %s. rank = %d\n", spasm_wtime() - start_time, nnz, u_n);
	
	free(qinv);
	free(p);
	return U;
}
