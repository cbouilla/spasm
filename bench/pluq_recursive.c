/* indent -nfbs -nip -npsl -di0 rank_hybrid.c */
#include <assert.h>
#include <stdio.h>
#include "spasm.h"
#include <getopt.h>

/** Computes the rank of the input matrix using the hybrid strategy */
int main(int argc, char **argv)
{
	char nnz[6];
	int allow_transpose = 1;	/* transpose ON by default */
	int n_times = 3;
	int prime = 42013;

	/* options descriptor */
	struct option longopts[7] = {
		{"max-recursion", required_argument, NULL, 'm'},
		{"no-transpose", no_argument, NULL, 'a'},
		{"modulus", required_argument, NULL, 'p'},
		{NULL, 0, NULL, 0}
	};

	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'm':
			n_times = atoi(optarg);
			break;
		case 'a':
			allow_transpose = 0;
			break;
		case 'p':
			prime = atoi(optarg);
			break;
		default:
			printf("Unknown option\n");
			exit(1);
		}
	}
	argc -= optind;
	argv += optind;


	spasm_triplet *T = spasm_load_sms(stdin, prime);
	if (allow_transpose && (T->n < T->m)) {
		fprintf(stderr, "[rank] transposing matrix : ");
		fflush(stderr);
		double start_time = spasm_wtime();
		spasm_triplet_transpose(T);
		fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
	}
	spasm *A = spasm_compress(T);
//	spasm *A_start = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;
	spasm_human_format(spasm_nnz(A), nnz);
	fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, nnz);

	/* allocate room for U */
	spasm *U = spasm_csr_alloc(n, m, spasm_nnz(A), A->prime, SPASM_WITH_NUMERICAL_VALUES);
	int *Up = U->p;
	int *Uj = U->j;
	int *Ux = U->x;
	int u_nnz = 0;
	int u_n = 0;
	double start_time = spasm_wtime();

	int *p = spasm_malloc(n * sizeof(int));
	int *qinv = spasm_malloc(m * sizeof(int));

	for (int iteration = 0; iteration < n_times; iteration++) {
		
		/* find new pivots */
		int npiv = spasm_find_pivots(A, p, qinv);
		spasm_make_pivots_unitary(A, p, npiv);
	
		/* copy them into U */
		int *Ap = A->p;
		int *Aj = A->j;
		int *Ax = A->x;	
		for (int i = 0; i < npiv; i++) {
			Up[u_n] = u_nnz;
			int I = p[i];
			/* not enough room in U ? realloc twice the size */
			if (u_nnz + m > U->nzmax) {
				spasm_csr_realloc(U, 2 * U->nzmax + m);
				Uj = U->j;
				Ux = U->x;
			}
			for (int px = Ap[I]; px < Ap[I + 1]; px++) {
				Uj[u_nnz] = Aj[px];
				Ux[u_nnz] = Ax[px];
				u_nnz++;
			}
			u_n++;
		}
		Up[u_n] = u_nnz;

		/* compute schur complement, update matrix */
		spasm *B = spasm_schur(A, p, qinv, npiv, -1, 0);
		spasm_csr_free(A);
		A = B;
	}

	/* final step : GPLU on the remainder */
	spasm_lu *LU = spasm_LU(A, SPASM_IDENTITY_PERMUTATION, SPASM_DISCARD_L);
	
	/* append last pivots to U */
	spasm_csr_free(A);
	A = LU->U;
	
	int *Ap = A->p;
	int *Aj = A->j;
	int *Ax = A->x;	
	for (int i = 0; i < A->n; i++) {
		/* not enough room in U ? realloc twice the size */
		if (u_nnz + m > U->nzmax) {
			spasm_csr_realloc(U, 2 * U->nzmax + m);
			Uj = U->j;
			Ux = U->x;
		}

		Up[u_n] = u_nnz;
		for (int px = Ap[i]; px < Ap[i + 1]; px++) {
			Uj[u_nnz] = Aj[px];
			Ux[u_nnz] = Ax[px];
			u_nnz++;
		}
		u_n++;
	}
	Up[u_n] = u_nnz;
	U->n = u_n;
	spasm_free_LU(LU);

	double end_time = spasm_wtime();
	spasm_human_format(spasm_nnz(U), nnz);
	fprintf(stderr, "done in %.3f s. NNZ(U) = %s. rank = %d\n", end_time - start_time, nnz, u_n);

#if 0
	/* check */
	for (int j = 0; j < m; j++)
		qinv[j] = -1;
	
	for (int i = 0; i < u_n; i++) {		
		for (int px = Up[i]; px < Up[i + 1]; px++)
			assert(qinv[Uj[px]] == -1);
		qinv[Uj[Up[i]]] = i;
	}

	int *xj = spasm_malloc(3*m * sizeof(int));
  	spasm_vector_zero(xj, 3*m);

  	int *x = spasm_malloc(m * sizeof(spasm_GFp));
  	int *y = spasm_malloc(m * sizeof(spasm_GFp));
  	int *z = spasm_malloc(n * sizeof(spasm_GFp));

  	assert(U->m == A_start->m);

  	for (int i = 0; i < n; i++) {
  		spasm_vector_zero(x, m);
  		spasm_vector_zero(y, m);
  		spasm_vector_zero(z, n);
  		
  		// in principle x.U == A_start[i]
  		int top = spasm_sparse_forward_solve(U, A_start, i, xj, x, qinv);
  		/*printf("top = %d, m=%d\n", top, m);
  		for (int px = top; px < m; px++)
  			printf("x[%d] = %d\n", xj[px], x[xj[px]]);
  		assert(x[Uj[Up[i]]] == 1);*/
  		for (int px = top; px < m; px++) {
  			int j = xj[px];
  			if (qinv[j] >= 0)
  				z[qinv[j]] = x[j];
  		}
  		// check it: y = x.U
  		spasm_gaxpy(U, z, y);
  		spasm_scatter(A_start->j, A_start->x, A_start->p[i], A_start->p[i + 1], A_start->prime - 1, y, A_start->prime);

  		for (int j = 0; j < m; j++)
  			if (y[j]) {
  				printf("i = %d\ny[%d] == %d\n", i, j, y[j]);
  				exit(1);
  			}
  	}
#endif

	spasm_save_csr(stdout, U);
	spasm_csr_free(U);
	//spasm_save_csr(A_start);
	free(p);
	free(qinv);
	return 0;
}
