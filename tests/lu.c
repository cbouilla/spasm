#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"

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
	parse_command_line_options(argc, argv);
	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;

	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);

	opts.complete = 1;
	struct spasm_lu *fact = spasm_echelonize(A, &opts);
	int r = fact->r;
	assert(fact->L != NULL);
	assert(r == 0 || fact->p != NULL);
	assert(fact->Ltmp == NULL);
	assert(r == fact->L->m);
	assert(r == fact->U->n);
	assert(fact->U->m == m);
	assert(fact->L->n == n);
	assert(spasm_factorization_verify(A, fact, 1337));
	assert(spasm_factorization_verify(A, fact, 21011984));

	bool *pivotal_row = spasm_malloc(n * sizeof(*pivotal_row));
	bool *pivotal_col = spasm_malloc(m * sizeof(*pivotal_row));
	for (int i = 0; i < n; i++)
		pivotal_row[i] = 0;
	for (int j = 0; j < m; j++)
		pivotal_col[j] = 0;
	for (int k = 0; k < r; k++) {
		int i = fact->p[k];
		pivotal_row[i] = 1;
	}
	for (int j = 0; j < m; j++) {
		int i = fact->qinv[j];
		if (i >= 0)
			pivotal_col[j] = 1;
	}

	#pragma omp parallel
	{

	spasm_ZZp *x = malloc(n * sizeof(*x));
	spasm_ZZp *y = malloc(m * sizeof(*y));
	spasm_ZZp *u = malloc(n * sizeof(*u));
	spasm_ZZp *v = malloc(m * sizeof(*v));

	/* check that A == L*U */
	#pragma omp for
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			x[j] = 0;
			u[j] = 0;
		}
		for (int j = 0; j < m; j++) {
			y[j] = 0;
			v[j] = 0;
		}
		// printf("\ri=%d / %d\n", i, n);
		// fflush(stdout);
		x[i] = 1;

		spasm_xApy(x, A, y);     // y <- x*A
		spasm_xApy(x, fact->L, u); // u <- x*L
		spasm_xApy(u, fact->U, v); // v <- (x*L)*U

		for (int j = 0; j < m; j++) {
			// printf("# x*A[%4d] = %8d VS x*L[%4d] = %8d VS x*LU[%4d] = %8d\n", j, y[j], j, u[j], j, v[j]);
			if (y[j] != v[j])
				printf("mismatch on row %d (pivotal=%d), column %d (pivotal=%d)\n", 
					i, pivotal_row[i], j, pivotal_col[j]);
			assert(y[j] == v[j]);
		}
		printf("OK on row %d (pivotal=%d)\n", i, pivotal_row[i]);

	}
	free(x);
	free(y);
	free(u);
	free(v);
	}
	assert(spasm_factorization_verify(A, fact, 42));
	printf("ok - L*U == A\n");
	spasm_csr_free(A);
	spasm_lu_free(fact);
	return 0;
}
