#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>

#include "spasm.h"

/** Computes the rank of the input matrix using the hybrid strategy described in [PASCO'17] */
int main(int argc, char **argv)
{
	int allow_transpose = 1;	/* transpose ON by default */
	int prime = 42013;

	/* parse command-line options */
	struct option longopts[7] = {
		{"no-transpose", no_argument, NULL, 'a'},
		{"modulus", required_argument, NULL, 'p'},
		{NULL, 0, NULL, 0}
	};
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
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

	/* load matrix from standard input */
	spasm_triplet *T = spasm_load_sms(stdin, prime);
	if (allow_transpose && (T->n < T->m)) {
		fprintf(stderr, "[rank] transposing matrix : ");
		fflush(stderr);
		double start = spasm_wtime();
		spasm_triplet_transpose(T);
		fprintf(stderr, "%.1f s\n", spasm_wtime() - start);
	}
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;
	char hnnz[8];
	spasm_human_format(spasm_nnz(A), hnnz);
	fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, hnnz);

	int *qinv = spasm_malloc(m * sizeof(int));
	double start_time = spasm_wtime();
	spasm *U = spasm_echelonize(A, qinv, NULL);   /* NULL = default options */
	double end_time = spasm_wtime();
	fprintf(stderr, "done in %.3f s rank = %d\n", end_time - start_time, U->n);
	spasm_csr_free(A);
	spasm_csr_free(U);
	free(qinv);
	return 0;
}
