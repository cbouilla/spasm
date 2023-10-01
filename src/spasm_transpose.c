#include <stdlib.h>

#include "spasm.h"

struct spasm_csr *spasm_transpose(const struct spasm_csr *C, int keep_values)
{
	int m = C->m;
	int n = C->n;
	const i64 *Cp = C->p;
	const int *Cj = C->j;
	const spasm_ZZp *Cx = C->x;
	i64 prime = spasm_get_prime(C);

	/* allocate result */
	struct spasm_csr *T = spasm_csr_alloc(m, n, spasm_nnz(C), prime, keep_values && (Cx != NULL));
	i64 *Tp = T->p;
	int *Tj = T->j;
	spasm_ZZp *Tx = T->x;

	/* get workspace */
	i64 *w = spasm_calloc(m, sizeof(*w));

	/* compute column counts */
	for (int i = 0; i < n; i++)
		for (i64 px = Cp[i]; px < Cp[i + 1]; px++) {
			int j = Cj[px];
			w[j] += 1;
		}

	/* compute column pointers (in both Cp and w) */
	i64 sum = 0;
	for (int j = 0; j < m; j++) {
		Tp[j] = sum;
		sum += w[j];
		w[j] = Tp[j];
	}
	Tp[m] = sum;

	/* dispatch entries */
	for (int i = 0; i < n; i++) {
		for (i64 px = Cp[i]; px < Cp[i + 1]; px++) {
			int j = Cj[px];
			i64 py = w[j];
			Tj[py] = i;
			if (Tx != NULL)
				Tx[py] = Cx[px];
			w[j] += 1;
		}
	}
	free(w);
	return T;
}
