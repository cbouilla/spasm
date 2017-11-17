#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/*
 * Finds pivots without performing arithmetic operations (using the
 * Faugère-Lachartre heuristic and some other ideas) permute the matrix so
 * that the pivots are on the top-left. The result is sent to the standard
 * output
 */

int main() {
	int principal = 220728;

	spasm_triplet *T = spasm_load_sms(stdin, 11);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	
	spasm *B = spasm_submatrix(A, principal, n, 0, principal, 1);
	
	spasm_save_csr(stdout, B);
	spasm_csr_free(B);
	spasm_csr_free(A);
	return 0;
}
