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

	/* compute LU factorization with L matrix */
	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	opts.L = 1;
	struct spasm_lu *fact = spasm_echelonize(A, &opts);
	assert(spasm_factorization_verify(A, fact, 42));

	bool *ok = spasm_malloc(n * sizeof(*ok));
	struct spasm_csr *X = spasm_gesv(fact, A, ok);
	assert(X->n == n);
	assert(X->m == n);
	for (int i = 0; i < n; i++)
		printf("ok[%d] = %d\n", i, ok[i]);

	/* check XA == A */
	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	spasm_ZZp *y = spasm_malloc(n * sizeof(*y));
	spasm_ZZp *z = spasm_malloc(m * sizeof(*y));
	spasm_ZZp *b = spasm_malloc(m * sizeof(*b));

	spasm_prng_ctx ctx;
	spasm_prng_seed_simple(prime, 0, 0, &ctx);
	for (int i = 0; i < n; i++) {
		x[i] = spasm_prng_ZZp(&ctx);
		y[i] = 0;
	}
	for (int j = 0; j < m; j++) {
		z[j] = 0;
		b[j] = 0;
	}
	
	spasm_xApy(x, A, b);
	spasm_xApy(x, X, y);
	spasm_xApy(y, A, z);

	for (int j = 0; j < m; j++)
		if (b[j] != z[j]) {
			printf("not ok - gesv solver [incorrect solution found]\n");
			exit(1);
		}
	printf("ok\n");
	return 0;
}
