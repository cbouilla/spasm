#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
	assert(argc == 2);
	int test = atoi(argv[1]);

	// load matrix
	spasm_triplet *T = spasm_load_sms(stdin, 32003);
	spasm *U = spasm_compress(T);
	spasm_triplet_free(T);
	int n = U->n;
	int m = U->m;

	assert(n<= m); // upper-trapezoidal

	// load RHS
	T = spasm_triplet_alloc(1, m, 10, 32003, true);
	spasm_add_entry(T, 0, 0, 1);
	spasm_add_entry(T, 0, m / 2, 2);
	spasm_add_entry(T, 0, m - 1, 3);
	spasm *B = spasm_compress(T);
	spasm_triplet_free(T);

	int *xi = malloc(3*m * sizeof(int));
	spasm_vector_zero(xi, 3*m);

	int *x = malloc(m * sizeof(spasm_GFp));
	int *y = malloc(m * sizeof(spasm_GFp));
	spasm_vector_zero(x, m);
	spasm_vector_zero(y, m);

	int *pinv = malloc(m * sizeof(int));
	for (int i = 0; i < n; i++)
		pinv[i] = i;
	for(int i = n; i < m; i++)
		pinv[i] = -1;

	spasm_sparse_forward_solve(U, B, 0, xi, x, pinv);

	spasm_gaxpy(U, x, y);
	for (int i = n; i < m; i++)
		y[i] = (y[i] + x[i]) % B->prime;
	spasm_scatter(B->j, B->x, B->p[0], B->p[1], B->prime - 1, y, B->prime);

	for (int i = 0; i < m; i++) {
		if (y[i] != 0) {
			printf("not ok %d - sparse triangular U-solve\n", test);
			exit(0);
		}
	}

	printf("ok %d - sparse triangular U-solve\n", test);

	spasm_csr_free(U);
	spasm_csr_free(B);
	free(xi);
	free(pinv);
	free(x);
	free(y);
	return 0;
}
