#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

int main(int argc, char **argv)
{
	struct spasm_triplet *T = spasm_triplet_load(stdin, 32003, NULL);
	struct spasm_csr *G = spasm_compress(T);
	spasm_triplet_free(T);
	int n = G->n;
	int m = G->m;
	assert(n >= m);
	assert(spasm_is_lower_triangular(G));

	spasm_ZZp *x = malloc(n * sizeof(spasm_ZZp));
	spasm_ZZp *b = malloc(m * sizeof(spasm_ZZp));
	spasm_ZZp *y = malloc(m * sizeof(spasm_ZZp));
	for (int i = 0; i < m; i++) {
		b[i] = spasm_ZZp_init(G->field, rand());
		y[i] = b[i];
	}

	spasm_dense_back_solve(G, y, x, SPASM_IDENTITY_PERMUTATION);

	for (int i = 0; i < m; i++)
		y[i] = 0;

	spasm_xApy(x, G, y);
	for (int i = 0; i < m; i++) {
		if (y[i] != b[i]) {
			printf("not ok - dense back-substitution triangular solver [incorrect solution found]\n");
			exit(1);
		}
	}
	printf("ok -  dense back-substitution triangular solver\n");
	spasm_csr_free(G);
	free(x);
	free(y);
	free(b);
	return 0;
}
