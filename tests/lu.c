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
	spasm_triplet *T = spasm_load_sms(stdin, prime);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;

	spasm_ZZp *x = malloc(n * sizeof(*x));
	spasm_ZZp *y = malloc(m * sizeof(*y));
	spasm_ZZp *u = malloc(n * sizeof(*u));
	spasm_ZZp *v = malloc(m * sizeof(*v));

	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	opts.L = 1;
	spasm_lu *fact = spasm_echelonize(A, &opts);
	assert(fact->L != NULL);
	assert(fact->U->n == 0 || fact->Lqinv != NULL);
	assert(fact->Ltmp == NULL);
	assert(fact->U->n == fact->L->m);
	assert(fact->U->m == m);
	assert(fact->L->n == n);

	/* check that A == L*U */
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			x[j] = 0;
			u[j] = 0;
		}
		for (int j = 0; j < m; j++) {
			y[j] = 0;
			v[j] = 0;
		}
		x[i] = 1;

		spasm_xApy(x, A, y);     // y <- x*A
		spasm_xApy(x, fact->L, u); // u <- x*L
		spasm_xApy(u, fact->U, v); // v <- (x*L)*U

		for (int j = 0; j < m; j++)
			if (y[j] != v[j]) {
				printf("not ok - L*U == A (col %d)\n", j);
				exit(1);
			}
	}
	printf("ok - L*U == A\n");
	spasm_csr_free(A);
	spasm_lu_free(fact);
	free(x);
	free(y);
	free(u);
	free(v);
	return 0;
}
