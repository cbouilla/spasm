#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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
	parse_command_line_options(argc, argv);
 	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
 	struct spasm_csr *A = spasm_compress(T);
 	spasm_triplet_free(T);
 	int m = A->m;
        if (m == 0) {
                printf("SKIP --- empty matrix / useless\n");
                exit(EXIT_SUCCESS);
        }

 	struct spasm_lu *fact = spasm_echelonize(A, NULL);
	struct spasm_csr *U = fact->U;
	int *qinv = fact->qinv;
	int r = U->n;
	if (r == 0) {
		printf("rank zero # SKIP\n");
		exit(EXIT_SUCCESS);
	}
 	int *pinv = spasm_malloc(r * sizeof(*pinv));
 	for (int j = 0; j < r; j++)
 		pinv[j] = -1;
 	for (int j = 0; j < m; j++)
 		if (qinv[j] >= 0)
 			pinv[qinv[j]] = j;
 	for (int j = 0; j < r; j++)
 		assert(pinv[j] != -1);

	struct spasm_csr *Ut = spasm_transpose(U, true);
	assert(Ut->n == m);
	assert(Ut->m == r);
  	spasm_csr_free(A);
  	spasm_csr_free(U);
	
	// make RHS
	T = spasm_triplet_alloc(1, r, 10, prime, true);
	spasm_add_entry(T, 0, 0, 1);
	// spasm_add_entry(T, 0, r / 2, 2);
	// spasm_add_entry(T, 0, r - 1, 3);
	struct spasm_csr *B = spasm_compress(T);
	spasm_triplet_free(T);

	for (int i = 0; i < r; i++)
		printf("# pinv[%d] = %d\n", i, pinv[i]);

	// int i = pinv[0];
	// const i64 *Utp = Ut->p;
	// const int *Utj = Ut->j;
	// printf("# U[%d] = ", i);
	// for (i64 px = Utp[i]; px < Utp[i + 1]; px++)
	// 	printf("%d ", Utj[px]);
	// printf("\n");

	// triangular solve
	int *xi = spasm_malloc(3*r * sizeof(*xi));
	for (int j = 0; j < 3*r; j++)
		xi[j] = 0;
	spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
	for (int i = 0; i < m; i++)
		x[i] = 0;
	int top = spasm_sparse_triangular_solve(Ut, B, 0, xi, x, pinv);

	for (int px = top; px < r; px++)
		printf("# x[%d] = %d\n", xi[px], x[xi[px]]);

	for (int i = 0; i < m; i++)
		printf("# xx[%d] = %d\n", i, x[i]);

	// check
	spasm_ZZp *xx = spasm_malloc(m * sizeof(*xx));
	spasm_ZZp *yy = spasm_malloc(r * sizeof(*yy));
	for (int j = 0; j < m; j++)
		xx[j] = 0;
	for (int j = 0; j < r; j++)
		yy[j] = 0;
	for (int px = top; px < r; px++) {
		int i = xi[px];
		if (x[i] == 0)
			continue;
		int j = pinv[i];
		if (j < 0) {
			yy[i] = x[i];
		}
		else {
			xx[j] = x[i];
		}
	}
	spasm_xApy(xx, Ut, yy);
	for (int i = 0; i < r; i++)
		printf("# y[%d] = %d\n", i, yy[i]);

        spasm_scatter(B, 0, -1, yy);
	for (int i = 0; i < r; i++)
		printf("# y[%d] = %d\n", i, yy[i]);

        for (int i = 0; i < r; i++)
                if (yy[i] != 0) {
                        printf("not ok - sparse triangular transpose(U)-solve (y[%d] == %d)\n", i, yy[i]);
                        exit(1);
                }
        printf("ok - sparse triangular transpose(U)-solve\n");
}
