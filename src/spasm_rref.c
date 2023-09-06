#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "spasm.h"

/* 
 * build the RREF from an echelonized matrix. 
 * This code is similar to src/spasm_schur.c ---> in bad need of factorization
 */
spasm * spasm_rref(const spasm *U, const int *Uqinv, int *Rqinv)
{
	int n = U->n;
	int m = U->m;
	char hnnz[8];
	spasm_human_format(spasm_nnz(U), hnnz);
	fprintf(stderr, "[rref] start. U is %d x %d (%s nnz)\n", n, m, hnnz);
	double start_time = spasm_wtime();
	spasm *R = spasm_csr_alloc(n, m, spasm_nnz(U), U->prime, SPASM_WITH_NUMERICAL_VALUES);
	i64 *Rp = R->p;
	int *Rj = R->j;
	int *Rx = R->x;
	int Rn = 0;         /* #rows in R */
	i64 nnz = 0;        /* entries in R */
	int writing = 0;
	const i64 *Up = U->p;
	const int *Uj = U->j;

	#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * m * sizeof(int));
		spasm_vector_zero(xj, 3 * m);
		int tid = spasm_get_thread_num();
		int *qinv_local = spasm_malloc(m * sizeof(int));
		for (int j = 0; j < m; j++)
			qinv_local[j] = Uqinv[j];
	
		#pragma omp for schedule(guided)
	  	for (int i = 0; i < n; i++) {
	  		int pivot = Uj[Up[i]];
	  		assert(qinv_local[pivot] == i);
	  		qinv_local[pivot] = -1;
	  		int top = spasm_sparse_triangular_solve(U, U, i, xj, x, qinv_local);
	  		
			/* ensure R has the "pivot first" property */
			for (int px = top + 1; px < m; px++)
				if (xj[px] == pivot) {
					xj[px] = xj[top];
					xj[top] = pivot;
					break;
				}
			assert(xj[top] == pivot);

	  		/* count the NZ in the new row */
	  		int row_nz = 0;
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if ((qinv_local[j] < 0) && (x[j] != 0))
					row_nz += 1;
			}

			int local_i;
			i64 local_nnz;
			#pragma omp critical(rref)
			{
				/* not enough room in R ? realloc twice the size */
				if (nnz + m > R->nzmax) {
					/* wait until other threads stop writing into R */
					#pragma omp flush(writing)
					while (writing > 0) {
						#pragma omp flush(writing)
					}
					spasm_csr_realloc(R, 2 * R->nzmax + m);
					Rj = R->j;
					Rx = R->x;
				}
				/* save row k */
				local_i = Rn;
				Rn += 1;
				local_nnz = nnz;
				nnz += row_nz;
				/* this thread will write into R */
				#pragma omp atomic update
				writing += 1;
			}

			/* write the new row in R */
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if (qinv_local[j] < 0 && x[j] != 0) {
					Rj[local_nnz] = xj[px];
					Rx[local_nnz] = x[j];
					local_nnz += 1;
				}
			}
			Rp[local_i + 1] = local_nnz;
			qinv_local[pivot] = i;

			/* we're done writing */
			#pragma omp atomic update
			writing--;

			if (tid == 0) {
				char hnnz[8];
				spasm_human_format(nnz, hnnz);
	  			fprintf(stderr, "\rRREF: %d/%d, |R| = %s    ", Rn, n, hnnz);
	  			fflush(stderr);
	  		}
		}
	  	free(x);
		free(xj);
		free(qinv_local);

		#pragma omp for
		for (int j = 0; j < m; j++)
			Rqinv[j] = -1;

		#pragma omp for
		for (int i = 0; i < n; i++) {
			int px = Rp[i];
			int j = Rj[px];
			Rqinv[j] = i;
		}
	}
	fprintf(stderr, "\n");

	spasm_human_format(spasm_nnz(R), hnnz);
	fprintf(stderr, "[rref] done in %.1fs. NNZ(R) = %s\n", spasm_wtime() - start_time, hnnz);
	return R;
}