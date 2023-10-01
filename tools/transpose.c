#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main() {
	struct spasm_triplet *T = spasm_triplet_load(stdin, 42013, NULL);

	/* TODO : this is really ugly */
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	struct spasm_csr *A_t = spasm_transpose(A, true);
	spasm_csr_free(A);

	spasm_csr_save(A_t, stdout);

	spasm_csr_free(A_t);
	return 0;
}
