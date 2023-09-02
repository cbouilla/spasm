#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "spasm.h"

/* build the RREF. This code is similar to src/spasm_schur.c ---> in bad need of factorization */
spasm * spasm_rref(const spasm *U, const int *Uqinv, int *Rqinv)
{
	int n = U->n;
	int m = U->m;
	char hnnz[8];
	spasm_human_format(spasm_nnz(U), hnnz);
	fprintf(stderr, "[rref] start. U is %d x %d (%s nnz)\n", n, m, hnnz);
	double start_time = spasm_wtime();
	spasm *R = spasm_csr_alloc(n, m, spasm_nnz(U), U->prime, SPASM_WITH_NUMERICAL_VALUES);
	int *Rp = R->p;
	int *Rj = R->j;
	int *Rx = R->x;
	int Rn = 0;         /* #rows in R */
	int nnz = 0;        /* entries in R */
	int writing = 0;
	const int *Up = U->p;
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
	  		int top = spasm_sparse_forward_solve(U, U, i, xj, x, qinv_local);
	  		
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

			int local_i, local_nnz;
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

/* given an echelonized matrix, return a basis of its right kernel */
spasm * spasm_kernel(const spasm *R, const int *qinv)
{
	assert(qinv != NULL);
	int n = R->n;
	int m = R->m;
	assert(n <= m);
	int prime = R->prime;
	spasm *Rt = spasm_transpose(R, SPASM_WITH_NUMERICAL_VALUES);
	const int *Rtp = Rt->p;
	const int *Rtj = Rt->j;
	const spasm_GFp *Rtx = Rt->x;

	int *p = spasm_malloc(n * sizeof(*p));    /* what row of Rt (column of R) contains the pivot on col i.*/
	const int *Rp = R->p;
	const int *Rj = R->j;
	for (int i = 0; i < n; i++) {
		int px = Rp[i];
		int j = Rj[px];
		p[i] = j;
	}
	spasm *K = spasm_csr_alloc(m - n, m, spasm_nnz(R) - n + m - n, prime, SPASM_WITH_NUMERICAL_VALUES);
	K->n = 0;
	int *Kp = K->p;
	int *Kj = K->j;
	spasm_GFp *Kx = K->x;
	int nnz = 0;      /* #entries in K */
	for (int j = 0; j < m; j++) {
		if (qinv[j] >= 0)
			continue;           /* skip pivotal columns of R */
		Kj[nnz] = j;
		Kx[nnz] = prime - 1;
		nnz += 1;
		for (int px = Rtp[j]; px < Rtp[j + 1]; px++) {
			int i = Rtj[px];
			Kj[nnz] = p[i];
			Kx[nnz] = Rtx[px];
			nnz += 1;
		}
		K->n += 1;
		Kp[K->n] = nnz;
	}
	assert(K->n == m - n);
	assert(nnz == spasm_nnz(R) - n + m - n);
	free(p);
	spasm_csr_free(Rt);
	return K;
}