#include <assert.h>
#include <stdlib.h>    // rand

#include "spasm.h"

/*
 * Permutations matrices are represented by vectors.
 * 
 * p[k] = i means that P[k,i] = 1
 */


/*
 * x <-- P.b (or, equivalently, x <-- b.(P^{-1}), for dense vectors x and b;
 * p=NULL denotes identity.
 * 
 * This means that x[k] <--- b[ p[k] ]
 */
void spasm_pvec(const int *p, const spasm_GFp *b, spasm_GFp *x, int n)
{
	assert(x != NULL);
	assert(b != NULL);
	for (int i = 0; i < n; i++) {
		int j = (p != SPASM_IDENTITY_PERMUTATION) ? p[i] : i;
		x[i] = b[j];
	}
}

/*
 * x <--- (P^{-1}).b (or x <--- b.P), for dense vectors x and b; p=NULL
 * denotes identity.
 * 
 * This means that x[ p[k] ] <--- b[ k ]
 * 
 * The function is given p, not p^{-1}.
 */
void spasm_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n)
{
	assert(x != NULL);
	assert(b != NULL);

	for (int i = 0; i < n; i++) {
		int j = (p != SPASM_IDENTITY_PERMUTATION) ? p[i] : i;
		x[j] = b[i];
	}
}

/* compute the inverse permutation */
int *spasm_pinv(int const *p, int n)
{
	/* p = NULL denotes identity */
	if (p == NULL) {
		return NULL;
	}
	/* allocate result */
	int *pinv = spasm_malloc(n * sizeof(*pinv));
	/* invert the permutation */
	for (int k = 0; k < n; k++)
		pinv[p[k]] = k;
	return pinv;
}


/*
 * C = P.A.Q^-1 where P and Q^-1 are permutations of 0..n-1 and 0..m-1
 * respectively.
 */
spasm *spasm_permute(const spasm *A, const int *p, const int *qinv, int values)
{
	/* check inputs */
	assert(A != NULL);
	int n = A->n;
	int m = A->m;
	const i64 *Ap = A->p;
	const int *Aj = A->j;
	const spasm_GFp *Ax = A->x;

	/* alloc result */
	spasm *C = spasm_csr_alloc(n, m, A->nzmax, A->prime, values && (Ax != NULL));
	i64 *Cp = C->p;
	int *Cj = C->j;
	spasm_GFp *Cx = C->x;
	i64 nnz = 0;

	for (int i = 0; i < n; i++) {
		/* row i of C is row p[i] of A (denoted by j) */
		Cp[i] = nnz;
		int j = (p != NULL) ? p[i] : i;
		for (i64 t = Ap[j]; t < Ap[j + 1]; t++) {
			/* col j of A is col qinv[j] of C */
			int j = Aj[t];
			int jj = (qinv != NULL) ? qinv[j] : j;
			Cj[nnz] = jj;
			if (Cx != NULL)
				Cx[nnz] = Ax[t];
			nnz += 1;
		}
	}
	/* finalize the last row of C */
	Cp[n] = nnz;
	return C;
}

int *spasm_random_permutation(int n) {
	int i, *p;

	p = spasm_malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		p[i] = i;
	}
	for (i = n - 1; i > 0; i--) {
		spasm_swap(p, i, rand() % i);
	}

	return p;
}

/* in-place permute x[a:b] using p. Destroys p */
void spasm_range_pvec(int *x, int a, int b, int *p) {
	int i;

	for (i = 0; i < b - a; i++) {
		p[i] = x[a + p[i]];
	}
	for (i = 0; i < b - a; i++) {
		x[a + i] = p[i];
	}
}
