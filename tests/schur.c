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
 
	int *p = spasm_malloc(n * sizeof(*p));
	int *qinv = spasm_malloc(m * sizeof(*qinv));
	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	
	struct spasm_csr *U = spasm_csr_alloc(n, m, spasm_nnz(A), prime, true);
	U->n = 0;
	for (int j = 0; j < m; j++)
		qinv[j] = -1;

	struct spasm_lu fact;
	fact.U = U;
	fact.qinv = qinv;
	fact.L = NULL;
	fact.Ltmp = NULL;
	fact.p = NULL;

	int npiv = spasm_pivots_extract_structural(A, NULL, &fact, p, &opts);
	struct spasm_csr *S = spasm_schur(A, p + npiv, n - npiv, &fact, -1, NULL, NULL, NULL);
	int Sn = S->n;
	i64 *Sp = S->p;
	int *Sj = S->j;

	/* checking that nothing remains under the pivots when we don't want it */
	for (int i = 0; i < Sn; i++) {
		for (i64 px = Sp[i]; px < Sp[i + 1]; px++) {
			int j = Sj[px];
			if (qinv[j] >= 0) {
				printf("not ok - coeff (%d, %d) is below a pivot\n", i, Sj[px]);
				exit(1);
			} 
		}
	}
	printf("ok - elimination coeffs are really absent\n");
	spasm_csr_free(S);

	/* phase II: check the elimination coeffs */
	/*
	int *p_out = spasm_malloc(Sn * sizeof(int));
	S = spasm_schur(A, p, npiv, -1, 1, p_out);
	Sp = S->p;
	Sj = S->j;
	spasm_ZZp *Sx = S->x;
	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));

	abort = 0;
	for (int k = 0; k < Sn; k++) {
		for (int i = 0; i < n; i++)
			x[i] = 0;
		for (int j = 0; j < m; j++)
			y[j] = 0;
		for (int px = Sp[k]; px < Sp[k + 1]; px++) {
			int j = Sj[px];
			if (qinv[j] >= 0)
				x[qinv[j]] = Sx[px];
			else
				y[j] = Sx[px];
		}

		assert(x[p_out[k]] == 0);		
		x[p_out[k]] = -1;
		
		spasm_gaxpy(A, x, y);

		
		for (int j = 0; j < m; j++)
			abort |= (y[j] != 0);
		if (abort) {
			printf("not ok %d - could not reconstruct row %d of S (=row %d of A)\n", test, k, p_out[k]);
			break;
		}
	}
	spasm_csr_free(S);
	if (!abort)
		printf("ok %d - reconstruction successfull\n", test);
	printf("not ok %d # TODO not implemented\n", test);
	*/

	spasm_csr_free(A);
	exit(EXIT_SUCCESS);
}