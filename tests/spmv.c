#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main()
{
	struct spasm_triplet *T = spasm_triplet_load(stdin, 257, NULL);
	struct spasm_csr *C = spasm_compress(T);
	spasm_triplet_free(T);

	int n = C->n;
	assert(n < spasm_get_prime(C));
	spasm_ZZp *x = malloc(n * sizeof(*x));
	spasm_ZZp *y = malloc(n * sizeof(*y));
	for (int i = 0; i < n; i++) {
		x[i] = spasm_ZZp_init(C->field, i + 1);
		y[i] = 0;
	}
	spasm_xApy(x, C, y);
	for (int i = 0; i < n; i++)
		printf("%d\n", y[i]);
	
	spasm_csr_free(C);
	free(x);
	free(y);
	return 0;
}
