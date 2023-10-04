#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <err.h>

#include "spasm.h"

/*
 * Generates a proof that A has the given rank.
 * This requires a full LU factorization.
 *
 * Implements the ideas given in
 * "A New Interactive Certificate for Matrix Rank", by Wayne Eberly.
 * Technical report, University of Calgary, 2015.
 * http://prism.ucalgary.ca/bitstream/1880/50543/1/2015-1078-11.pdf. 
 *
 * To make a "real" certificate, seed should be obtained by hashing the input matrix
 * and the "commitment" i / j.
 */
struct spasm_rank_certificate * spasm_certificate_rank_create(const struct spasm_csr *A, const u8 *hash, const struct spasm_lu *fact)
{
	assert(fact->L != NULL);
	const struct spasm_csr *U = fact->U;
	const struct spasm_csr *L = fact->L;
	int n = L->n;
	int m = U->m;
	int r = U->n;

	struct spasm_rank_certificate *proof = spasm_malloc(sizeof(*proof));
	proof->r = r;
	for (int i = 0; i < 32; i++)
		proof->hash[i] = hash[i];
	proof->prime = spasm_get_prime(A);
	spasm_ZZp *xx = spasm_malloc(r * sizeof(*xx));
	spasm_ZZp *yy = spasm_malloc(r * sizeof(*yy));
	int *ii = spasm_malloc(r * sizeof(*ii));
	int *jj = spasm_malloc(r * sizeof(*jj));
	proof->x = xx;
	proof->y = yy;
	proof->i = ii;
	proof->j = jj;

	/* write i / j indices (positions of pivots) */
	for (int k = 0; k < r; k++)
		ii[k] = fact->p[k];
	int k = 0;
	for (int j = 0; j < m; j++)
		if (fact->qinv[j] >= 0) {
			jj[k] = j;
			k += 1;
		}

	/* Generate challenge */
	spasm_prng_ctx ctx;
	spasm_prng_seed(hash, proof->prime, 0, &ctx);
	
	/* compute x */
	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	for (int j = 0; j < m; j++)
		y[j] = 0;
	for (int k = 0; k < r; k++) {
		int j = jj[k];
		y[j] = spasm_prng_ZZp(&ctx);
	}
	spasm_solve(fact, y, x);
	for (int k = 0; k < r; k++) {
		int i = ii[k];
		xx[k] = x[i];
	}
	
	/* compute y */
	spasm_ZZp BOT = 0x7fffffff;
	for (int i = 0; i < n; i++)
		x[i] = BOT;
	for (int k = 0; k < r; k++) {
		int i = ii[k];
		x[i] = 0;
	}
	for (int i = 0; i < n; i++)
		if (x[i] == BOT)
			x[i] = -spasm_prng_ZZp(&ctx);
	for (int j = 0; j < m; j++)
		y[j] = 0;
	spasm_xApy(x, A, y);
	spasm_solve(fact, y, x);
	for (int k = 0; k < r; k++) {
		int i = ii[k];
		yy[k] = x[i];
	}
	free(x);
	free(y);
	return proof;
}

bool spasm_certificate_rank_verify(const struct spasm_csr *A, const u8 *hash, const struct spasm_rank_certificate *proof)
{
	int n = A->n;
	int m = A->m;
	int r = proof->r;

	for (int i = 0; i < 32; i++)
		if (hash[i] != proof->hash[i])
			return 0;

	if (spasm_get_prime(A) != proof->prime)
		return 0;

	for (int k = 0; k < r; k++) {
		if (proof->i[k] < 0 || proof->i[k] >= n)
			return 0;
		if (proof->j[k] < 0 || proof->j[k] >= m)
			return 0;
	}

	spasm_prng_ctx ctx;
	spasm_prng_seed(proof->hash, proof->prime, 0, &ctx);
	spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	
	/* check Ax --- must match alpha on the j coordinates */
	for (int i = 0; i < n; i++)   /* compute x */
		x[i] = 0;
	for (int k = 0; k < r; k++) {
		int i = proof->i[k];
		x[i] = proof->x[k];
	}
	for (int j = 0; j < m; j++)
		y[j] = 0;
	spasm_xApy(x, A, y);
	bool correct = 1;
	for (int k = 0; k < r; k++) {
		int j = proof->j[k];
		if (y[j] != spasm_prng_ZZp(&ctx))
			correct = 0;
	}

	/* check Ay */
	spasm_ZZp BOT = 0x7fffffff;        /* compute y */
	for (int i = 0; i < n; i++)
		x[i] = BOT;
	for (int k = 0; k < r; k++) {
		int i = proof->i[k];
		x[i] = proof->y[k];
	}
	//int k = 0;
	for (int i = 0; i < n; i++)
		if (x[i] == BOT)
			x[i] = spasm_prng_ZZp(&ctx);
	for (int j = 0; j < m; j++)
		y[j] = 0;
	spasm_xApy(x, A, y);
	for (int j = 0; j < m; j++) {
		if (y[j] != 0)
			correct = 0;
	}
	free(x);
	free(y);
	return correct;
}

/* TODO: move this out of here */
bool spasm_factorization_verify(const struct spasm_csr *A, const struct spasm_lu *fact, u64 seed)
{
	assert(fact->L != NULL);
	const struct spasm_csr *U = fact->U;
	const struct spasm_csr *L = fact->L;
	// const int *Uqinv = fact->qinv;
	const int *Lp = fact->p;

	int n = A->n;
	int m = A->m;
	int r = U->n;
	spasm_ZZp *x = malloc(n * sizeof(*x));
	spasm_ZZp *y = malloc(r * sizeof(*y));
	spasm_ZZp *z = malloc(m * sizeof(*z));
	spasm_ZZp *t = malloc(m * sizeof(*t));
	bool correct = 1;
	bool complete = 0;
	
	bool *pivotal_row = spasm_malloc(n * sizeof(*pivotal_row));
	for (int i = 0; i < n; i++)
		pivotal_row[i] = 0;
	for (int j = 0; j < r; j++) {
		int i = Lp[j];
		assert(i >= 0);
		pivotal_row[i] = 1;
	}

	spasm_prng_ctx ctx;
	spasm_prng_seed_simple(spasm_get_prime(A), seed, 0, &ctx);
	for (int j = 0; j < m; j++) {
		z[j] = 0;
		t[j] = 0;
	}
	for (int k = 0; k < r; k++)
		y[k] = 0;
	for (int i = 0; i < n; i++) {
		spasm_ZZp foo = spasm_prng_ZZp(&ctx);
		if (complete || pivotal_row[i])
			x[i] = foo;
		else
			x[i] = 0;
	}
	spasm_xApy(x, A, t);
	spasm_xApy(x, L, y);
	spasm_xApy(y, U, z);
	for (int j = 0; j < m; j++)
		if (z[j] != t[j])
			correct = 0;

	free(pivotal_row);
	free(x);
	free(y);
	free(z);
	free(t);
	return correct;
}

void spasm_rank_certificate_save(const struct spasm_rank_certificate *proof, FILE *f)
{
	int r = proof->r;
	fprintf(f, "%d\n", r);
	fprintf(f, "%" PRId64 "\n", proof->prime);
	for (int i = 0; i < 32; i++)
		fprintf(f, "%02x", proof->hash[i]);
	fprintf(f, "\n");
	for (int k = 0; k < r; k++)
		fprintf(f, "%d ", proof->i[k]);
	fprintf(f, "\n");
	for (int k = 0; k < r; k++)
		fprintf(f, "%d ", proof->j[k]);
	fprintf(f, "\n");
	for (int k = 0; k < r; k++)
		fprintf(f, "%d ", proof->x[k]);
	fprintf(f, "\n");
	for (int k = 0; k < r; k++)
		fprintf(f, "%d ", proof->y[k]);
	fprintf(f, "\n");
}

bool spasm_rank_certificate_load(FILE *f, struct spasm_rank_certificate *proof)
{
	int r;
	if (1 != fscanf(f, "%d", &r))
		return 0;
	proof->r = r;
	proof->i = malloc(r * sizeof(*proof->i));
	proof->j = malloc(r * sizeof(*proof->i));
	proof->x = malloc(r * sizeof(*proof->i));
	proof->y = malloc(r * sizeof(*proof->i));
	if (1 != fscanf(f, "%" SCNd64 "\n", &proof->prime))
		return 0;
	char hash[65];
	if (NULL == fgets(hash, 65, f))
		return 0;
	for (int i = 0; i < 32; i++) {
		char byte[3] = {hash[2*i], hash[2*i+1], 0};
		proof->hash[i] = strtoul(byte, NULL, 16);
	}
	for (int k = 0; k < r; k++)
		fscanf(f, "%d", &proof->i[k]);
	for (int k = 0; k < r; k++)
		fscanf(f, "%d", &proof->i[k]);
	for (int k = 0; k < r; k++)
		fscanf(f, "%d", &proof->x[k]);
	for (int k = 0; k < r; k++)
		fscanf(f, "%d", &proof->y[k]);
	return 1;
}