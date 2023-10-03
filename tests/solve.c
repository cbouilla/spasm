#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"
#include "test_tools.h"

i64 prime = 42013;

void parse_command_line_options(int argc, char **argv)
{
	struct option longopts[] = {
		{"modulus", required_argument, NULL, 'p'},
		{NULL, 0, NULL, 0}
	};
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'p':
			prime = atoll(optarg);
			break;
		default:
			errx(1, "Unknown option\n");
		}
	}
}

int main(int argc, char **argv)
{
	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;

	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	spasm_ZZp *b = spasm_malloc(m * sizeof(*b));

	/* compute LU factorization with L matrix */
	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	opts.L = 1;
	struct spasm_lu *fact = spasm_echelonize(A, &opts);
	assert(spasm_factorization_verify(A, fact, 42));

	/* test solve() with a sensible RHS */
	printf("# testing correct solution\n");
	
	/* forge valid solution */
	spasm_prng_ctx ctx;
	spasm_prng_seed_simple(prime, 0, 0, &ctx);
	for (int i = 0; i < n; i++)
		x[i] = spasm_prng_ZZp(&ctx);
	for (int j = 0; j < m; j++)
		b[j] = 0;
	spasm_xApy(x, A, b);

	for (int i = 0; i < n; i++) /* forget solution */
		x[i] = 0;

	bool result = spasm_solve(fact, b, x);  /* find solution */
	if (!result) {
		printf("not ok - LU solver [solution not found]\n");
		exit(1);
	}

	for (int j = 0; j < m; j++)   /* check solution */
		y[j] = 0;
	spasm_xApy(x, A, y);
	for (int j = 0; j < m; j++)
		if (y[j] != b[j]) {
			printf("not ok - LU solver [incorrect solution found]\n");
			exit(1);
		}

	/* test with a bogus RHS */
	if (n < m) {
		printf("# testing bogus solution\n");
		for (int j = 0; j < m; j++)       /* create random vector */
			b[j] = spasm_prng_ZZp(&ctx);
	
		bool result = spasm_solve(fact, b, x);
		if (result) {
			printf("not ok - LU solver [bogus solution found]\n");
			exit(1);
		}

	}
	printf("ok - LU solver\n");
	spasm_csr_free(A);
	spasm_lu_free(fact);
	free(x);
	free(y);
	free(b);
	return 0;
}
