#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

/* 
 * given the transpose of an echelonized matrix (not necessarily in RREF), return a basis of its right kernel
 * This code is similar to src/spasm_schur.c and src/spasm_rref.c ---> in bad need of factorization
 * qinv locates the pivots in Ut.
 */
spasm * spasm_kernel(const spasm *Ut, const int *qinv)
{
	int m = Ut->n;
	int n = Ut->m;
	int prime = Ut->prime;
	assert(n <= m);
	char hnnz[8];
	spasm_human_format(spasm_nnz(Ut), hnnz);
	fprintf(stderr, "[kernel] start. U is %d x %d (%s nnz)\n", n, m, hnnz);
	double start_time = spasm_wtime();
	spasm *K = spasm_csr_alloc(m-n, m, spasm_nnz(Ut), prime, SPASM_WITH_NUMERICAL_VALUES);
	i64 *Kp = K->p;
	int *Kj = K->j;
	spasm_GFp *Kx = K->x;
	int Kn = 0;         /* #rows in R */
	i64 nnz = 0;        /* entries in R */
	int writing = 0;

	int *Utqinv = spasm_malloc(n * sizeof(*Utqinv));
	for (int j = 0; j < m; j++) {
		int i = qinv[j];
		if (i >= 0) {
			assert(i < n);
			Utqinv[i] = j;
		}
	}

	#pragma omp parallel
	{
		spasm_GFp *x = spasm_malloc(m * sizeof(*x));
		spasm_GFp *y = spasm_malloc(n * sizeof(*x));
		int *xj = spasm_malloc(3 * n * sizeof(int));
		spasm_vector_zero(xj, 3 * n);
		int tid = spasm_get_thread_num();

		#pragma omp for schedule(guided)
	  	for (int j = 0; j < m; j++) {
	  		if (qinv[j] >= 0)
	  			continue;         /* skip pivotal row */
	  		int top = spasm_sparse_backward_solve(Ut, Ut, j, xj, x, Utqinv);

	  		/* debug */
	  		for (int k = top; k < n; k++) {
	  			int j = xj[k];
	  			assert(j < m);
	  			assert(qinv[j] >= 0);
	  		}
	  		/* debug: check the product */
			spasm_vector_zero(y, n);
			spasm_gaxpy(Ut, x, y);
        	spasm_scatter(Ut->j, Ut->x, Ut->p[j], Ut->p[j + 1], prime - 1, y, prime);
        	for (int i = 0; i < n; i++)
				assert(y[i] == 0);

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
					fprintf(stderr, "increasing K; %" PRId64 " ---> %" PRId64 "\n", K->nzmax, 2 * K->nzmax + row_nz);
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
			Kx[local_nnz] = prime - 1;
			local_nnz += 1;
			for (int px = top; px < n; px++) {
				int j = xj[px];
				if (x[j] != 0) {
					Kj[local_nnz] = xj[px];
					Kx[local_nnz] = x[j];
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
	  			fprintf(stderr, "\rkernel: %d/%d, |K| = %s    ", Kn, m-n, hnnz);
	  			fflush(stderr);
	  		}
		}
	  	free(x);
		free(xj);
	}
	fprintf(stderr, "\n");
	free(Utqinv);
	spasm_human_format(spasm_nnz(K), hnnz);
	fprintf(stderr, "[kernel] done in %.1fs. NNZ(K) = %s\n", spasm_wtime() - start_time, hnnz);
	return K;
}


/* given an echelonized matrix (in RREF), return a basis of its right kernel */
spasm * spasm_kernel_from_rref(const spasm *R, const int *qinv)
{
	assert(qinv != NULL);
	int n = R->n;
	int m = R->m;
	assert(n <= m);
	int prime = R->prime;
	spasm *Rt = spasm_transpose(R, SPASM_WITH_NUMERICAL_VALUES);
	const i64 *Rtp = Rt->p;
	const int *Rtj = Rt->j;
	const spasm_GFp *Rtx = Rt->x;

	int *p = spasm_malloc(n * sizeof(*p));    /* what row of Rt (column of R) contains the pivot on col i.*/
	const i64 *Rp = R->p;
	const int *Rj = R->j;
	for (int i = 0; i < n; i++) {
		int px = Rp[i];
		int j = Rj[px];
		p[i] = j;
	}
	spasm *K = spasm_csr_alloc(m - n, m, spasm_nnz(R) - n + m - n, prime, SPASM_WITH_NUMERICAL_VALUES);
	K->n = 0;
	i64 *Kp = K->p;
	int *Kj = K->j;
	spasm_GFp *Kx = K->x;
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