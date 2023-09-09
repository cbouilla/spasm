#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv) 
{
	int prime = 42013;

	spasm_triplet *T = spasm_load_sms(stdin, prime);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	int n = A->n;
	int m = A->m;
  
	int *p = spasm_malloc(n * sizeof(int));
	int *qinv = spasm_malloc(m * sizeof(int));
	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	int npiv = spasm_find_pivots(A, p, qinv, &opts);
	spasm_make_pivots_unitary(A, p, npiv);

	spasm *S = spasm_schur(A, p, npiv, A, qinv, -1, SPASM_DISCARD_L, NULL);
	int Sn = S->n;
	i64 *Sp = S->p;
	int *Sj = S->j;

	/* checking that nothing remains under the pivots when we don't want it */
	for (int i = 0; i < Sn; i++) {
		for (i64 px = Sp[i]; px < Sp[i + 1]; px++)
			if (qinv[Sj[px]] >= 0) {
				printf("not ok - coeff (%d, %d) is below a pivot\n", i, Sj[px]);
				exit(1);
			} 
	}
	printf("ok - elimination coeffs are really absent\n");
	spasm_csr_free(S);

	/* phase II: check the elimination coeffs */
	/*
	int *p_out = spasm_malloc(Sn * sizeof(int));
	S = spasm_schur(A, p, npiv, -1, 1, p_out);
	Sp = S->p;
	Sj = S->j;
	spasm_GFp *Sx = S->x;
	spasm_GFp *x = spasm_malloc(n * sizeof(*x));
	spasm_GFp *y = spasm_malloc(m * sizeof(*y));

	abort = 0;
	for (int k = 0; k < Sn; k++) {
		for (int i = 0; i < n; i++)
			x[i] = 0;
		for (int j = 0; j < m; j++)
			y[j] = 0;
		for (int px = Sp[k]; px < Sp[k + 1]; px++) {
			int j = Sj[px];
			if (qinv[j] >= 0)
				x[qinv[j]] = Sx[px];
			else
				y[j] = Sx[px];
		}

		assert(x[p_out[k]] == 0);		
		x[p_out[k]] = prime - 1;
		
		spasm_gaxpy(A, x, y);

		
		for (int j = 0; j < m; j++)
			abort |= (y[j] != 0);
		if (abort) {
			printf("not ok %d - could not reconstruct row %d of S (=row %d of A)\n", test, k, p_out[k]);
			break;
		}
	}
	spasm_csr_free(S);
	if (!abort)
		printf("ok %d - reconstruction successfull\n", test);
	printf("not ok %d # TODO not implemented\n", test);
	*/

	spasm_csr_free(A);
	exit(EXIT_SUCCESS);
}