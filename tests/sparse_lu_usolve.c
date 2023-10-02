#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"
#include "test_tools.h"
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
 	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
 	struct spasm_csr *A = spasm_compress(T);
 	spasm_triplet_free(T);
 	int m = A->m;
        if (m <= 0) {
                printf("SKIP --- empty matrix / useless\n");
                exit(EXIT_SUCCESS);
        }

 	struct spasm_lu *fact = spasm_echelonize(A, NULL);
 	struct spasm_csr *U = fact->U;
 	int *qinv = fact->qinv;
	int r = U->n;
 	
	// make RHS
	T = spasm_triplet_alloc(1, m, 10, prime, true);
	spasm_add_entry(T, 0, 0, 1);
	if (m > 0) {
		spasm_add_entry(T, 0, m / 2, 2);
		spasm_add_entry(T, 0, m - 1, 3);
	}
	struct spasm_csr *B = spasm_compress(T);
	spasm_triplet_free(T);

	// triangular solve
	int *xj = malloc(3*m * sizeof(*xj));
	for (int j = 0; j < 3*m; j++)
                xj[j] = 0;
	spasm_ZZp *x = malloc(m * sizeof(spasm_ZZp));
	int top = spasm_sparse_triangular_solve(U, B, 0, xj, x, qinv);

	for (int px = top; px < m; px++) {
		int j = xj[px];
		printf("# x[%d] = %d / qinv[%d] = %d\n", j, x[j], j, qinv[j]);
	}

	spasm_ZZp *xx = malloc(r * sizeof(*xx));
	spasm_ZZp *yy = malloc(m * sizeof(*yy));
	for (int j = 0; j < r; j++)
                xx[j] = 0;
        for (int j = 0; j < m; j++)
                yy[j] = 0;
	for (int px = top; px < m; px++) {
		int j = xj[px];
		if (x[j] == 0)
			continue;
		int i = qinv[j];
		if (i < 0) {
			yy[j] = x[j];
			printf("y[%d] <-- %d\n", j, x[j]);
		}
		else {
			xx[i] = x[j];
			printf("x[%d] <-- %d\n", i, x[j]);
		}
	}
	spasm_xApy(xx, U, yy);
	for (int i = 0; i < m; i++)
		printf("# y[%d] = %d\n", i, yy[i]);

        spasm_scatter(B, 0, -1, yy);
	for (int i = 0; i < m; i++)
		printf("# y[%d] = %d\n", i, yy[i]);

        for (int i = 0; i < m; i++)
                if (yy[i] != 0) {
                        printf("not ok - sparse LU / triangular-solve (y[%d] == %d)\n", i, yy[i]);
                        exit(1);
                }
        printf("ok - sparse triangular LU / U-solve\n");
}
