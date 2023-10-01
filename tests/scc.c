#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv)
{
	struct spasm_triplet *T = spasm_triplet_load(stdin, 42013, NULL);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);

	int n = A->n;
	int m = A->m;
	if (n != m) {
		printf("ok - # SKIP because matrix is non-square");
		exit(0);
	}

	/* generate random row & col permutation */
	/*p = spasm_random_permutation(n);
	pinv = spasm_pinv(p, n);
	B = spasm_permute(A, p, pinv, true);
	free(pinv);
	spasm_csr_free(A);
	*/
	struct spasm_csr *B = A;
	struct spasm_dm *P = spasm_strongly_connected_components(B);
	int *p = P->p;
	int *q = P->q;
	int *rr = P->r;
	int nb = P->nb;

	/* verbosity */
	printf("# p = ");
	for (int k = 0; k < n; k++)
		printf("%d ", p[k]);
 
	printf("\n# -----------------------\n");
	printf("#rr = ");
	for (int k = 0; k < nb + 1; k++)
		printf("%d ", rr[k]);
	printf("\n ");

	/* --- check that p and q match ---------- */
	for(int i = 0; i < n; i++)
	       if (p[i] != q[i]) {
			printf("not ok - SCC - p != q\n");
			exit(1);
		}
  
	/* --- check that p is actually a permutation ---------------- */
	int *x = spasm_malloc(n * sizeof(*x));
	for (int i = 0; i < n; i++)
		x[i] = 0;
	for(int i = 0; i < n; i++) {
		if (p[i] < 0 || p[i] >= n) {
			printf("not ok - SCC - p is out of range, p[%d] = %d\n", i, p[i]);
			exit(1);
		}
		x[ p[i] ]++;
	}
	for (int i = 0; i < n; i++)
		if (x[i] != 1) {
			printf("not ok - SCC - p is not bijective (x[%d] = %d)\n", i, x[i]);
			exit(1);
		}
	free(x);

	/* --- verbosity ---------------- */
	printf("# SCC = %d\n", nb);

	for (int k = 0; k < nb; k++) {
		printf("# SCC_%d : ", k);
		for (int i = rr[k]; i < rr[k + 1]; i++)
			printf("%d ", p[i]);
		printf("\n");
	}

	/* --- check that decomposition is really block-upper-triangular ---------------- */
	int *pinv = spasm_pinv(p, n);
	struct spasm_csr *C = spasm_permute(B, p, pinv, SPASM_IGNORE_VALUES);
	i64 *Cp = C->p;
	int *Cj = C->j;
	free(pinv);
	spasm_csr_free(B);

	for (int k = 0; k < nb; k++)
		for (int i = rr[k]; i < rr[k + 1]; i++)
			for (i64 px = Cp[i]; px < Cp[i + 1]; px++) {
				int j = Cj[px];
				if (j < rr[k]) {
					printf("not ok - SCC - row %d (in C_%d) has entries on column %d\n", i, k, j);
					exit(1);
				}
			}
 	printf("ok - SCC\n");
	spasm_csr_free(C);
	return 0;
}
