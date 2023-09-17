#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <err.h>

#include "spasm.h"

/*
 * Implements the ideads given in 
 * "Elimination-based certificates for triangular equivalence and rank profiles"
 * By Jean-Guillaume Dumas, Erich Kaltofen, David Lucas and Clément Pernet
 * J. Symb. Comp. --- https://doi.org/10.1016/j.jsc.2019.07.013
 */


/*
 * Generates a proof that U is correct.
 * This requires a full LU factorization.
 * To make a "real" certificate, seed should be obtained by hashing the input matrix
 */
spasm_rowspan_certificate * spasm_certificate_rowspan_create(const spasm *A, const spasm_lu *fact, u64 seed)
{
	assert(fact->L != NULL);
	const spasm *U = fact->U;
	const spasm *L = fact->L;
	int n = L->n;
	int m = U->m;
	int r = U->n;

	spasm_rowspan_certificate *proof = spasm_malloc(sizeof(*proof));
	proof->seed = seed;
	spasm_prng_seed(seed, 0);

	/* 
	 * Generate random vector x, then compute z s.t. xA == zU
	 * z "proves" that rowspan(A) is contained into rowspan(U)
	 */
	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	spasm_ZZp *z = spasm_malloc(r * sizeof(*z));
	for (int i = 0; i < n; i++)
		x[i] = spasm_ZZp_init(A->field, spasm_prng_next());
	for (int j = 0; j < m; j++)
		y[j] = 0;
	spasm_xApy(x, A, y);

	/* inverse permutation for U */
	int *Uq = spasm_malloc(r * sizeof(*Uq));
	const int *Uqinv = fact->Uqinv;
	for (int j = 0; j < m; j++) {
		int i = Uqinv[j];
		if (i != -1)
			Uq[i] = j;
	}
	bool ok = spasm_dense_forward_solve(U, y, z, Uq);
	assert(ok);
	free(x);
	free(y);
	free(Uq);
	proof->z = z;

	/*
	 * Generate random vector f, then compute t s.t. fU == hA
	 * h "proves" that rowspan(U) is contained into rowspan(A)
	 */
	spasm_ZZp *f = spasm_malloc(r * sizeof(*x));
	spasm_ZZp *g = spasm_malloc(m * sizeof(*y));
	spasm_ZZp *h = spasm_malloc(n * sizeof(*z));
	for (int i = 0; i < r; i++)
		f[i] = spasm_ZZp_init(A->field, spasm_prng_next());
	for (int j = 0; j < m; j++)
		g[j] = 0;
	spasm_xApy(f, U, g);
	ok = spasm_solve(fact, g, h);
	assert(ok);
	free(f);
	free(g);
	proof->h = h;
	return proof;
}

/* check that U is actually permuted upper-triangular */
static bool check_echelon_form(const spasm *U)
{
	int m = U->m;
	int r = U->n;
	const i64 *Up = U->p;
	const int *Uj = U->j;
	const spasm_ZZp *Ux = U->x;

	int *w = spasm_malloc(m * sizeof(*w));
	int *v = spasm_malloc(r * sizeof(*v));
	for (int j = 0; j < m; j++)
		w[j] = -1;
	/* mark columns with the last row where they occur */
	for (int i = 0; i < r; i++)
		for (i64 px = Up[i]; px < Up[i + 1]; px++) {
			if (Ux[px] == 0)
				continue;
			int j = Uj[px];
			w[j] = i;
		}
	/* bucket sort */
	for (int i = 0; i < r; i++)
		v[i] = 0;
	for (int j = 0; j < m; j++) {
		int i = w[j];
		if (i >= 0)
			v[i] = 1;
	}
	/* check that each row appear on a column of its own */
	bool correct = 1;
	for (int i = 0; i < r; i++)
		if (v[i] == 0)
			correct = 0;
	free(w);
	free(v);
	return correct;
}

/* check that xA == zB */
static bool check_vector_equality(const spasm_ZZp *x, const spasm *A, const spasm_ZZp *z, const spasm *B)
{
	int m = A->m;
	assert(m == B->m);
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	for (int j = 0; j < m; j++)
		y[j] = 0;
	spasm_xApy(x, A, y);
	for (int j = 0; j < m; j++)
		y[j] = -y[j];
	spasm_xApy(z, B, y);
	bool correct = 1;
	for (int j = 0; j < m; j++)
		if (y[j] != 0)
			correct = 0;
	free(y);
	return correct;
}

/*
 * Generates a proof that U is correct.
 * This requires a full LU factorization.
 *
 * The proof only requires U to be checked, and the seed.
 * To make a "real" certificate, seed should be obtained by hashing the input matrix
 */
bool spasm_certificate_rowspan_verify(const spasm *A, const spasm *U, const spasm_rowspan_certificate *proof)
{
	if (!check_echelon_form(U))
		return 0;

	int n = A->n;
	int r = U->n;
	
	/* check that xA == zU */
	spasm_prng_seed(proof->seed, 0);
	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	for (int i = 0; i < n; i++)
		x[i] = spasm_ZZp_init(A->field, spasm_prng_next());
	bool correct = check_vector_equality(x, A, proof->z, U);

	/* check that xU == hA */
	for (int i = 0; i < r; i++)
		x[i] = spasm_ZZp_init(A->field, spasm_prng_next());
	correct = correct && check_vector_equality(x, U, proof->h, A);
	free(x);
	return correct;
}