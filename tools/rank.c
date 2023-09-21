#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>

#include "spasm.h"

i64 prime = 42013;
bool compute_cert = 0;
bool allow_transpose = 1;	/* transpose ON by default */
struct echelonize_opts opts;

void parse_command_line_options(int argc, char **argv)
{
	struct option longopts[7] = {
		{"no-transpose", no_argument, NULL, 'a'},
		{"modulus", required_argument, NULL, 'p'},
		{"certificate", no_argument, NULL, 'c'},
		{NULL, 0, NULL, 0}
	};
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'a':
			allow_transpose = 0;
			break;
		case 'p':
			prime = atoll(optarg);
			break;
		case 'c':
			compute_cert = 1;
			break;
		default:
			printf("Unknown option\n");
			exit(1);
		}
	}
}


/** Computes the rank of the input matrix using the hybrid strategy described in [PASCO'17] */
int main(int argc, char **argv)
{
	spasm_echelonize_init_opts(&opts);
	parse_command_line_options(argc, argv);

	/* load matrix from standard input */
	spasm_triplet *T = spasm_load_sms(stdin, prime);
	if (allow_transpose && (T->n < T->m)) {
		fprintf(stderr, "[rank] transposing matrix : ");
		fflush(stderr);
		spasm_triplet_transpose(T);
	}
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;
	char hnnz[8];
	spasm_human_format(spasm_nnz(A), hnnz);
	fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, hnnz);

	if (compute_cert)
		opts.L = 1;

	double start_time = spasm_wtime();
	spasm_lu *fact = spasm_echelonize(A, &opts);   /* NULL = default options */
	double end_time = spasm_wtime();
	fprintf(stderr, "done in %.3f s rank = %d\n", end_time - start_time, fact->U->n);
	

	spasm_rank_certificate *proof = NULL;
	if (compute_cert) {
		assert(spasm_factorization_verify(A, fact, 42));
		assert(spasm_factorization_verify(A, fact, 1337));
		assert(spasm_factorization_verify(A, fact, 21011984));
		u64 seed = 42;
		fprintf(stderr, "generating certificate\n");
		proof = spasm_certificate_rank_create(A, fact, seed);
		fprintf(stderr, "checking certificate\n");
		bool correct = spasm_certificate_rank_verify(A, proof);
		if (correct)
			fprintf(stderr, "CORRECT certificate\n");
		else
			fprintf(stderr, "INCORRECT certificate\n");
	}
	spasm_lu_free(fact);
	spasm_csr_free(A);
	return 0;
}
