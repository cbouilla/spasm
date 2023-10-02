#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <err.h>
#include <inttypes.h>

#include "spasm.h"

int spasm_get_num_threads() {
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

int spasm_get_thread_num() {
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}


double spasm_wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}


i64 spasm_nnz(const struct spasm_csr * A)
{
	return A->p[A->n];
}

/* return a string representing n in 8 bytes */
void spasm_human_format(i64 n, char *target)
{
	if (n < 1000) {
		sprintf(target, "%" PRId64, n);
		return;
	}
	if (n < 1000000) {
		sprintf(target, "%.1fk", n / 1e3);
		return;
	}
	if (n < 1000000000) {
		sprintf(target, "%.1fm", n / 1e6);
		return;
	}
	if (n < 1000000000000ll) {
		sprintf(target, "%.1fg", n / 1e9);
		return;
	}
	if (n < 1000000000000000ll) {
		sprintf(target, "%.1ft", n / 1e12);
		return;
	}
}

void *spasm_malloc(i64 size)
{
	void *x = malloc(size);
	if (x == NULL)
		err(1, "malloc failed (size %" PRId64 ")", size);
	return x;
}

void *spasm_calloc(i64 count, i64 size)
{
	void *x = calloc(count, size);
	if (x == NULL)
		err(1, "calloc failed");
	return x;
}

void *spasm_realloc(void *ptr, i64 size)
{
	void *x = realloc(ptr, size);
	if (ptr != NULL && x == NULL && size != 0)
		err(1, "realloc failed");
	return x;
}

/* allocate a sparse matrix (compressed-row form) */
struct spasm_csr *spasm_csr_alloc(int n, int m, i64 nzmax, i64 prime, bool with_values)
{
	struct spasm_csr *A = spasm_malloc(sizeof(*A));
	spasm_field_init(prime, A->field);
	A->m = m;
	A->n = n;
	A->nzmax = nzmax;
	A->p = spasm_malloc((n + 1) * sizeof(i64));
	A->j = spasm_malloc(nzmax * sizeof(int));
	A->x = with_values ? spasm_malloc(nzmax * sizeof(spasm_ZZp)) : NULL;
	A->p[0] = 0;
	return A;
}

/* allocate a sparse matrix (triplet form) */
struct spasm_triplet *spasm_triplet_alloc(int n, int m, i64 nzmax, i64 prime, bool with_values)
{
	struct spasm_triplet *A = spasm_malloc(sizeof(*A));
	A->m = m;
	A->n = n;
	A->nzmax = nzmax;
	spasm_field_init(prime, A->field);
	A->nz = 0;
	A->i = spasm_malloc(nzmax * sizeof(int));
	A->j = spasm_malloc(nzmax * sizeof(int));
	A->x = with_values ? spasm_malloc(nzmax * sizeof(spasm_ZZp)) : NULL;
	return A;
}

/*
 * change the max # of entries in a sparse matrix. If nzmax < 0, then the
 * matrix is trimmed to its current nnz.
 */
void spasm_csr_realloc(struct spasm_csr *A, i64 nzmax)
{
	if (nzmax < 0)
		nzmax = spasm_nnz(A);
	// if (spasm_nnz(A) > nzmax)
	// 	errx(1, "spasm_csr_realloc with too small nzmax (contains %" PRId64 " nz, asking nzmax=%" PRId64 ")", spasm_nnz(A), nzmax);
	A->j = spasm_realloc(A->j, nzmax * sizeof(int));
	if (A->x != NULL)
		A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_ZZp));
	A->nzmax = nzmax;
}

/*
 * change the max # of entries in a sparse matrix. If nzmax < 0, then the
 * matrix is trimmed to its current nnz.
 */
void spasm_triplet_realloc(struct spasm_triplet *A, i64 nzmax)
{
	if (nzmax < 0)
		nzmax = A->nz;
	// fprintf(stderr, "[realloc] nzmax=%ld. before %px %px %px\n", nzmax, A->i, A->j, A->x);
	A->i = spasm_realloc(A->i, nzmax * sizeof(int));
	A->j = spasm_realloc(A->j, nzmax * sizeof(int));
	if (A->x != NULL)
		A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_ZZp));
	// fprintf(stderr, "[realloc] after %px %px %px\n", A->i, A->j, A->x);
	A->nzmax = nzmax;
}

/* free a sparse matrix */
void spasm_csr_free(struct spasm_csr *A)
{
	if (A == NULL)
		return;
	free(A->p);
	free(A->j);
	free(A->x);		/* trick : free does nothing on NULL pointer */
	free(A);
}

void spasm_triplet_free(struct spasm_triplet *A)
{
	free(A->i);
	free(A->j);
	free(A->x);		/* trick : free does nothing on NULL pointer */
	free(A);
}

void spasm_csr_resize(struct spasm_csr *A, int n, int m)
{
	A->m = m;
	/* TODO: in case of a shrink, check that no entries are left outside */
	A->p = spasm_realloc(A->p, (n + 1) * sizeof(i64));
	if (A->n < n) {
		i64 *Ap = A->p;
		for (int i = A->n; i < n + 1; i++)
			Ap[i] = Ap[A->n];
	}
	A->n = n;
}

struct spasm_dm * spasm_dm_alloc(int n, int m)
{
	struct spasm_dm *P = spasm_malloc(sizeof(*P));
	P->p = spasm_malloc(n * sizeof(int));
	P->q = spasm_malloc(m * sizeof(int));
	P->r = spasm_malloc((n + 6) * sizeof(int));
	P->c = spasm_malloc((m + 6) * sizeof(int));
	P->nb = 0;
	for (int i = 0; i < 5; i++) {
		P->rr[i] = 0;
		P->cc[i] = 0;
	}
	return P;
}

void spasm_dm_free(struct spasm_dm *P)
{
	free(P->p);
	free(P->q);
	free(P->r);
	free(P->c);
	free(P);
}

void spasm_lu_free(struct spasm_lu *N)
{
	free(N->qinv);
	free(N->p);
	spasm_csr_free(N->U);
	spasm_csr_free(N->L);
	free(N);
}
