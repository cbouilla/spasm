#include <stdio.h>
#include <stdlib.h>

#include "spasm.h"

/*
 * Finds pivots without performing arithmetic operations (using the
 * Faugère-Lachartre heuristic and some other ideas) permute the matrix so
 * that the pivots are on the top-left. The result is sent to the standard
 * output
 */

int main() {
	int prime = 42013;
	spasm_triplet *T = spasm_load_sms(stdin, prime);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;

	int *p = spasm_malloc(n * sizeof(int));
	int *qinv = spasm_malloc(m * sizeof(int));
	int npiv = spasm_find_pivots(A, p, qinv, NULL);
	spasm *B = spasm_permute_pivots(A, p, qinv, npiv);

	spasm_save_csr(stdout, B);
	free(p);
	free(qinv);
	spasm_csr_free(B);
	spasm_csr_free(A);
	return 0;
}
