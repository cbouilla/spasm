#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv)
{
 	int prime = 32003;
 	// load matrix
 	spasm_triplet *T = spasm_load_sms(stdin, prime);
 	spasm *A = spasm_compress(T);
 	spasm_triplet_free(T);
 	int m = A->m;

 	int *qinv = spasm_malloc(m * sizeof(*qinv));
 	spasm *U = spasm_echelonize(A, qinv, NULL);
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

	spasm *Ut = spasm_transpose(U, SPASM_WITH_NUMERICAL_VALUES);
	assert(Ut->n == m);
	assert(Ut->m == r);
  	spasm_csr_free(A);
  	spasm_csr_free(U);
	
	// make RHS
	T = spasm_triplet_alloc(1, r, 10, prime, true);
	spasm_add_entry(T, 0, 0, 1);
	// spasm_add_entry(T, 0, r / 2, 2);
	// spasm_add_entry(T, 0, r - 1, 3);
	spasm *B = spasm_compress(T);
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
	int *xi = malloc(3*r * sizeof(*xi));
	for (int j = 0; j < 3*r; j++)
		xi[j] = 0;
	spasm_GFp *x = malloc(m * sizeof(*x));
	int top = spasm_sparse_triangular_solve(Ut, B, 0, xi, x, pinv);

	for (int px = top; px < r; px++)
		printf("# x[%d] = %d\n", xi[px], x[xi[px]]);

	for (int i = 0; i < m; i++)
		printf("# xx[%d] = %d\n", i, x[i]);

	// check
	spasm_GFp *xx = malloc(m * sizeof(spasm_GFp));
	spasm_GFp *yy = malloc(r * sizeof(spasm_GFp));
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

        spasm_scatter(B, 0, prime - 1, yy);
	for (int i = 0; i < r; i++)
		printf("# y[%d] = %d\n", i, yy[i]);

        for (int i = 0; i < r; i++)
                if (yy[i] != 0) {
                        printf("not ok - sparse triangular transpose(U)-solve (y[%d] == %d)\n", i, yy[i]);
                        exit(1);
                }
        printf("ok - sparse triangular transpose(U)-solve\n");
}
