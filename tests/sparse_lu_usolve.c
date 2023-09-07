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
 	
	// make RHS
	T = spasm_triplet_alloc(1, m, 10, prime, true);
	spasm_add_entry(T, 0, 0, 1);
	if (m > 0) {
		spasm_add_entry(T, 0, m / 2, 2);
		spasm_add_entry(T, 0, m - 1, 3);
	}
	spasm *B = spasm_compress(T);
	spasm_triplet_free(T);

	// triangular solve
	int *xj = malloc(3*m * sizeof(*xj));
	spasm_vector_zero(xj, 3*m);
	spasm_GFp *x = malloc(m * sizeof(spasm_GFp));
	int top = spasm_sparse_triangular_solve(U, B, 0, xj, x, qinv);

	for (int px = top; px < m; px++) {
		int j = xj[px];
		printf("# x[%d] = %d / qinv[%d] = %d\n", j, x[j], j, qinv[j]);
	}

	spasm_GFp *xx = malloc(r * sizeof(spasm_GFp));
	spasm_GFp *yy = malloc(m * sizeof(spasm_GFp));
	spasm_vector_zero(xx, r);
	spasm_vector_zero(yy, m);
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

        spasm_scatter(B->j, B->x, B->p[0], B->p[1], prime - 1, yy, prime);
	for (int i = 0; i < m; i++)
		printf("# y[%d] = %d\n", i, yy[i]);

        for (int i = 0; i < m; i++)
                if (yy[i] != 0) {
                        printf("not ok - sparse LU / triangular-solve (y[%d] == %d)\n", i, yy[i]);
                        exit(1);
                }
        printf("ok - sparse triangular LU / U-solve\n");
}
