#include <assert.h>
#include <stdlib.h>

#include "spasm.h"

/* add an entry to a triplet matrix; enlarge it if necessary */
void spasm_add_entry(spasm_triplet *T, int i, int j, spasm_GFp x)
{
	assert((i >= 0) && (j >= 0));
	int prime = T->prime;

	i64 px = T->nz;
	if (px == T->nzmax)
		spasm_triplet_realloc(T, 1 + 2 * T->nzmax);
	if (T->x != NULL) {
		spasm_GFp x_p = ((x % prime) + prime) % prime;
		if (x_p == 0)
			return;
		T->x[px] = x_p;
	}

	T->i[px] = i;
	T->j[px] = j;
	T->nz += 1;
	T->n = spasm_max(T->n, i + 1);
	T->m = spasm_max(T->m, j + 1);
}

void spasm_triplet_transpose(spasm_triplet *T)
{
	int *foo = T->i;
	T->i = T->j;
	T->j = foo;
	int bar = T->m;
	T->m = T->n;
	T->n = bar;
}

/* in-place */
void spasm_deduplicate(spasm *A)
{
	int m = A->m;
	int n = A->n;
	i64 *Ap = A->p;
	int *Aj = A->j;
	spasm_GFp *Ax = A->x;
	int prime = A->prime;

	i64 *v = spasm_malloc(m * sizeof(*v));
	for (int j = 0; j < m; j++)
		v[j] = -1;

	i64 nz = 0;
	for (int i = 0; i < n; i++) {
		i64 p = nz;
		for (i64 it = Ap[i]; it < Ap[i + 1]; it++) {
			int j = Aj[it];
			if (v[j] < p) { /* occurs in previous row */
				v[j] = nz;
				Aj[nz] = j;
				if (Ax)
					Ax[nz] = Ax[it];
				nz += 1;
			} else {
				if (Ax)
					Ax[v[j]] = (Ax[v[j]] + Ax[it]) % prime;
			}
		}
		Ap[i] = p;
	}
	Ap[n] = nz;
	free(v);
	spasm_csr_realloc(A, -1);
}

/* C = compressed-row form of a triplet matrix T */
spasm *spasm_compress(const spasm_triplet * T)
{
	int m = T->m;
	int n = T->n;
	i64 nz = T->nz;
	int *Ti = T->i;
	int *Tj = T->j;
	spasm_GFp *Tx = T->x;
	
	double start = spasm_wtime();
	fprintf(stderr, "[CSR] Compressing... ");
	fflush(stderr);

	/* allocate result */
	spasm *C = spasm_csr_alloc(n, m, nz, T->prime, Tx != NULL);

	/* get workspace */
	i64 *w = spasm_calloc(n, sizeof(*w));
	i64 *Cp = C->p;
	int *Cj = C->j;
	spasm_GFp *Cx = C->x;

	/* compute row counts */
	for (int it = 0; it < nz; it++) {
		int i = Ti[it];
		w[i] += 1;
	}

	/* compute row pointers (in both Cp and w) */
	i64 sum = 0;
	for (int k = 0; k < n; k++) {
		Cp[k] = sum;
		sum += w[k];
		w[k] = Cp[k];
	}
	Cp[n] = sum;

	/* dispatch entries */
	for (i64 k = 0; k < nz; k++) {
		int i = Ti[k];
		i64 px = w[i];
		w[i] += 1;
		Cj[px] = Tj[k];
		if (Cx != NULL)
			Cx[px] = Tx[k];
	}
	free(w);
	spasm_deduplicate(C);

	/* success; free w and return C */
	char mem[16];
	int size = sizeof(int) * (n + nz) + sizeof(spasm_GFp) * ((Cx != NULL) ? nz : 0);
	spasm_human_format(size, mem);
	fprintf(stderr, "%" PRId64 " actual NZ, Mem usage = %sbyte [%.2fs]\n", spasm_nnz(C), mem, spasm_wtime() - start);
	return C;
}