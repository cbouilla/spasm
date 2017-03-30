/* indent -nfbs -nip -npsl -di0 matching_pivots.c */
#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/*
 * Finds pivots without performing arithmetic operations (using the
 * Faugère-Lachartre heuristic and some other ideas) and computes the Schur
 * complement w.r.t. these pivots. The result is sent to the standard output
 */

int main() {
	int npiv, n, m, prime = 42013, *p, *qinv;
	spasm_triplet *T;
	spasm *A;
	//double schur_density;

	T = spasm_load_sms(stdin, prime);
	A = spasm_compress(T);
	spasm_triplet_free(T);
	n = A->n;
	m = A->m;

	p = spasm_malloc(n * sizeof(int));
	qinv = spasm_malloc(m * sizeof(int));

	npiv = spasm_pivots_ur_matching(A, p, qinv);
	A = spasm_permute_pivots(A, p, qinv, npiv);
	//spasm_make_pivots_unitary(A, p, npiv);

	/* /\* estimate an upper-bound on the rank of the complement *\/ */
	/* schur_density = spasm_schur_probe_density(A, p, qinv, npiv, 100); */
	/* int nnz = schur_density * (n - npiv) * (m - npiv); */
	/* char tmp[6]; */
	/* spasm_human_format(sizeof(int) * (n - npiv + nnz) + sizeof(spasm_GFp) * nnz, tmp); */
	/* fprintf(stderr, "Schur complement: (%d x %d), estimted density : %.4f (%s byte)\n", n - npiv, m - npiv, schur_density, tmp); */

	/* /\* go for it *\/ */
	/* S = spasm_schur(A, p, qinv, npiv); */
	/* fprintf(stderr, "Schur complement: (%d x %d), nnz : %d, dens : %.5f\n", S->n, S->m, spasm_nnz(S), 1. * spasm_nnz(S) / (1. * S->n * S->m)); */

	 spasm_save_csr(stdout, A); 
	free(p);
	free(qinv);
	//spasm_csr_free(S);
	spasm_csr_free(A);
	return 0;
}
