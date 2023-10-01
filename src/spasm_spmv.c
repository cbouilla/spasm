#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

/*
 * (dense) vector * (sparse) Matrix
 * y <--- x*A + y
 */
void spasm_xApy(const spasm_ZZp *x, const struct spasm_csr *A, spasm_ZZp *y)
{
	int n = A->n;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_ZZp *Ax = A->x;
	for (int i = 0; i < n; i++)
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			y[j] = spasm_ZZp_axpy(A->field, x[i], Ax[px], y[j]);
		}
}

/*
 * (sparse) Matrix * (dense) vector
 * y <--- A*x + y
 */
void spasm_Axpy(const struct spasm_csr *A, const spasm_ZZp *x, spasm_ZZp *y)
{
	int n = A->n;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_ZZp *Ax = A->x;
	for (int i = 0; i < n; i++)
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			y[i] = spasm_ZZp_axpy(A->field, x[j], Ax[px], y[i]);
		}
}
