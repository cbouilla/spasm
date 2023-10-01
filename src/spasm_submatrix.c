#include <assert.h>
#include "spasm.h"

/**
 * returns A[r_0:r_1, c_0:c_1]
 */
struct spasm_csr *spasm_submatrix(const struct spasm_csr * A, int r_0, int r_1, int c_0, int c_1, int with_values)
{
	assert(A != NULL);
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_ZZp *Ax = A->x;
	i64 prime = spasm_get_prime(A);

	int Bn = spasm_max(0, r_1 - r_0);
	int Bm = spasm_max(0, c_1 - c_0);
	i64 Bnz = spasm_max(0, Ap[r_1] - Ap[r_0]);
	struct spasm_csr *B = spasm_csr_alloc(Bn, Bm, Bnz, prime, (A->x != NULL) && with_values);
	i64 *Bp = B->p;
	int *Bj = B->j;
	spasm_ZZp *Bx = B->x;

	i64 k = 0;
	for (int i = r_0; i < r_1; i++) {
		Bp[i - r_0] = k;
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			if (c_0 <= j && j < c_1) {
				Bj[k] = j - c_0;
				if (Bx != NULL)
					Bx[k] = Ax[px];
				k += 1;
			}
		}
	}

	/* finalize */
	Bp[r_1 - r_0] = k;

	/* shrink */
	spasm_csr_realloc(B, -1);
	return B;
}
