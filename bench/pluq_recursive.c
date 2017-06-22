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
	spasm *A_start = spasm_compress(T);
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

	spasm_human_format(spasm_nnz(U), nnz);
	fprintf(stderr, "LU factorization done in %.3f s. NNZ(U) = %s. rank = %d\n", spasm_wtime() - start_time, nnz, u_n);

	start_time = spasm_wtime();

	/* check that U is permuted lower-triangular*/
	for (int j = 0; j < m; j++)
		qinv[j] = -1;
	
	for (int i = 0; i < u_n; i++) {		
		for (int px = Up[i]; px < Up[i + 1]; px++)
			assert(qinv[Uj[px]] == -1);
		qinv[Uj[Up[i]]] = i;
	}

	/* build p to order the rows */
	int k = 0;
	for (int j = 0; j < m; j++)
		if (qinv[j] >= 0)
			p[k++] = qinv[j];
	assert(k == u_n);

	/* build the RREF */
	spasm *R = spasm_csr_alloc(u_n, m, spasm_nnz(U), A->prime, SPASM_WITH_NUMERICAL_VALUES);
	int *Rp = R->p;
	int *Rj = R->j;
	int *Rx = R->x;	
	int rnz = 0;

	int *xj = spasm_malloc(3*m * sizeof(int));
  	int *x = spasm_malloc(m * sizeof(spasm_GFp));
  	spasm_vector_zero(xj, 3*m);

  	for (int i = 0; i < u_n; i++) {
  		spasm_human_format(rnz, nnz);
  		fprintf(stderr, "\rRREF: %d/%d, |R| = %s    ", i, u_n, nnz);
  		fflush(stderr);

  		int j = Uj[Up[i]];
  		assert(qinv[j] == i);
  		qinv[j] = -1;
  		int top = spasm_sparse_forward_solve(U, U, i, xj, x, qinv);
  		
		/* not enough room in R ? realloc twice the size */
		if (rnz + m > R->nzmax) {
			spasm_csr_realloc(R, 2 * R->nzmax + m);
			Rj = R->j;
			Rx = R->x;
		}

		/* ensure R has the "pivot first" property */
		for (int px = top + 1; px < m; px++)
			if (xj[px] == j) {
				xj[px] = xj[top];
				xj[top] = j;
				break;
			}
		assert(xj[top] == j);

		/* copy into R */
		Rp[i] = rnz;
		for (int px = top; px < m; px++) {
			int j = xj[px];
			if (qinv[j] < 0) {
				Rj[rnz] = xj[px];
				Rx[rnz] = x[xj[px]];
				rnz++;
			}
		}

	}
	Rp[u_n] = rnz;
	fprintf(stderr, "\n");
	

	/* check */
	for (int j = 0; j < m; j++)
		qinv[j] = -1;
	for (int i = 0; i < u_n; i++)
		qinv[Uj[Up[i]]] = i;

	for (int i = 0; i < u_n; i++) {
		assert(qinv[Rj[Rp[i]]] == i);
		int found = 0;
		for (int px = Rp[i]; px < Rp[i + 1]; px++) {
			found |= Rj[px] == Uj[Up[i]];
			if (qinv[Rj[px]] != -1)
				assert(qinv[Rj[px]] == i);
		}
		assert(found);
	}
	spasm_csr_free(U);

	for (int i = 0; i < u_n; i++)
		qinv[Rj[Rp[p[i]]]] = i;

	spasm *S = spasm_permute(R, p, SPASM_IDENTITY_PERMUTATION, SPASM_WITH_NUMERICAL_VALUES);
	free(p);
	spasm_csr_free(R);

	spasm_human_format(spasm_nnz(S), nnz);
	fprintf(stderr, "done in %.3f s. NNZ(R) = %s\n", spasm_wtime() - start_time, nnz);
	
	spasm_save_csr(stdout, S);
	spasm_csr_free(S);
	//spasm_save_csr(A_start);
	free(qinv);
	free(x);
	free(xj);
	return 0;
}
