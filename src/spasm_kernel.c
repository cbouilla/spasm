#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

/* 
 * return a basis of the right kernel of the matrix described by the LU factorization
 */
struct spasm_csr * spasm_kernel(const struct spasm_lu *fact)
{
	const struct spasm_csr *U = fact->U;
	const int *qinv = fact->qinv; 
	int m = U->m;
	int n = U->n;
	i64 prime = spasm_get_prime(U);
	assert(n <= m);
	char hnnz[8];
	spasm_human_format(spasm_nnz(U), hnnz);
	logprintf("[kernel] start. U is %d x %d (%s nnz). Transposing U\n", n, m, hnnz);
	double start_time = spasm_wtime();

	struct spasm_csr *Ut = spasm_transpose(U, true);
	
	struct spasm_csr *K = spasm_csr_alloc(m-n, m, spasm_nnz(U), prime, true);
	i64 *Kp = K->p;
	int *Kj = K->j;
	spasm_ZZp *Kx = K->x;
	int Kn = 0;         /* #rows in R */
	i64 nnz = 0;        /* entries in R */
	int writing = 0;

	int *Utqinv = spasm_malloc(n * sizeof(*Utqinv)); /* locate pivots in Ut */
	for (int j = 0; j < m; j++) {
		int i = qinv[j];
		if (i >= 0)
			Utqinv[i] = j;
	}

	/*
	 * The following code is simiar to spasm_schur and spasm_rref (needs factorization...)
	 */
	#pragma omp parallel
	{
		spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
		int *xj = spasm_malloc(3 * n * sizeof(int));
		for (int j = 0; j < 3 * n; j++)
			xj[j] = 0;
		int tid = spasm_get_thread_num();

		#pragma omp for schedule(guided)
	  	for (int j = 0; j < m; j++) {
	  		if (qinv[j] >= 0)
	  			continue;         /* skip pivotal row */
	  		int top = spasm_sparse_triangular_solve(Ut, Ut, j, xj, x, Utqinv);

	  		/* count the NZ in the new row */
	  		int row_nz = 1;
	  		for (int px = top; px < n; px++) {
				int j = xj[px];
				if (x[j] != 0)
					row_nz += 1;
			}

	  		/* commit the row */
			int local_i;
			i64 local_nnz;
			#pragma omp critical(kernel)
			{
				/* not enough room in R ? realloc twice the size */
				if (nnz + row_nz > K->nzmax) {
					/* wait until other threads stop writing into R */
					#pragma omp flush(writing)
					while (writing > 0) {
						#pragma omp flush(writing)
					}
					logprintf("increasing K; %" PRId64 " ---> %" PRId64 "\n", K->nzmax, 2 * K->nzmax + row_nz);
					spasm_csr_realloc(K, 2 * K->nzmax + row_nz);
					Kj = K->j;
					Kx = K->x;
				}
				/* save row k */
				local_i = Kn;
				Kn += 1;
				local_nnz = nnz;
				nnz += row_nz;
				/* this thread will write into K */
				#pragma omp atomic update
				writing += 1;
			}

			/* write the new row in K */
			Kj[local_nnz] = j;
			Kx[local_nnz] = -1;
			local_nnz += 1;
			for (int px = top; px < n; px++) {
				int jj = xj[px];
				if (x[jj] != 0) {
					int ii = Utqinv[jj];
					Kj[local_nnz] = ii;
					Kx[local_nnz] = x[jj];
					local_nnz += 1;
				}
			}
			Kp[local_i + 1] = local_nnz;

			/* we're done writing */
			#pragma omp atomic update
			writing--;

			if (tid == 0) {
				char hnnz[8];
				spasm_human_format(nnz, hnnz);
	  			logprintf("\rkernel: %d/%d, |K| = %s    ", Kn, m-n, hnnz);
	  			fflush(stderr);
	  		}
		}
	  	free(x);
		free(xj);
	}

	logprintf("\n");
	free(Utqinv);
	spasm_csr_free(Ut);
	spasm_human_format(spasm_nnz(K), hnnz);
	logprintf("[kernel] done in %.1fs. NNZ(K) = %s\n", spasm_wtime() - start_time, hnnz);
	return K;
}

/* 
 * given an echelonized matrix (in RREF), return a basis of its right kernel.
 * This is less computationnaly expensive than from a non-RREF matrix
 */
struct spasm_csr * spasm_kernel_from_rref(const struct spasm_csr *R, const int *qinv)
{
	assert(qinv != NULL);
	int n = R->n;
	int m = R->m;
	assert(n <= m);
	i64 prime = spasm_get_prime(R);
	struct spasm_csr *Rt = spasm_transpose(R, true);
	const i64 *Rtp = Rt->p;
	const int *Rtj = Rt->j;
	const spasm_ZZp *Rtx = Rt->x;

	int *p = spasm_malloc(n * sizeof(*p));    /* what row of Rt (column of R) contains the pivot on col i.*/
	const i64 *Rp = R->p;
	const int *Rj = R->j;
	for (int i = 0; i < n; i++) {
		int px = Rp[i];
		int j = Rj[px];
		p[i] = j;
	}
	struct spasm_csr *K = spasm_csr_alloc(m - n, m, spasm_nnz(R) - n + m - n, prime, true);
	K->n = 0;
	i64 *Kp = K->p;
	int *Kj = K->j;
	spasm_ZZp *Kx = K->x;
	i64 nnz = 0;      /* #entries in K */
	for (int j = 0; j < m; j++) {
		if (qinv[j] >= 0)
			continue;           /* skip pivotal columns of R */
		Kj[nnz] = j;
		Kx[nnz] = prime - 1;
		nnz += 1;
		for (i64 px = Rtp[j]; px < Rtp[j + 1]; px++) {
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
