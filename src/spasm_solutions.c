#include <assert.h>
#include <stdlib.h>

#include "spasm.h"

/*
 * Solves x.A = b where A is not necessarily square;
 * 
 * b has size m, solution has size n.
 * 
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 * 
 * Uses a PLUQ factorization
 */
int spasm_PLUQ_solve(spasm * A, const spasm_ZZp * b, spasm_ZZp * x) {
	spasm_ZZp *u, *v, *w, *s;
	spasm_lu *PLUQ;
	spasm *L, *U;
	int n, m, r, ok;
	int *p, *qinv;

	/* check inputs */
	assert(A != NULL);
	assert(b != NULL);

	n = A->n;
	m = A->m;
	p = spasm_malloc(n * sizeof(int));
	qinv = spasm_malloc(m * sizeof(int));

	spasm_find_pivots(A, p, qinv);
	PLUQ = spasm_PLUQ(A, p, SPASM_KEEP_L);
	L = PLUQ->L;
	U = PLUQ->U;
	r = U->n;

	/* get workspace */
	u = spasm_malloc(m * sizeof(spasm_ZZp));
	v = spasm_malloc(r * sizeof(spasm_ZZp));
	w = spasm_malloc(n * sizeof(spasm_ZZp));
	s = spasm_malloc(n * sizeof(spasm_ZZp));

	/* u*Q = b */
	spasm_ipvec(PLUQ->qinv, b, u, m);

	/* v.U*Q = b  (if possible) */
	ok = spasm_dense_forward_solve(U, u, v, SPASM_IDENTITY_PERMUTATION);

	if (ok == SPASM_SUCCESS) {
		/* w.LUQ = b */
		spasm_dense_back_solve(L, v, w, SPASM_IDENTITY_PERMUTATION);

		/* x.PLUQ = b */
		spasm_ipvec(PLUQ->p, w, s, n);
		spasm_ipvec(p, s, x, n);
	}
	free(u);
	free(v);
	free(w);
	free(s);
	free(p);
	free(qinv);
	spasm_free_LU(PLUQ);
	return ok;
}

/*
 * Solves x.A = b where A is not necessarily square;
 * 
 * b has size m, solution has size n.
 * 
 * returns SPASM_SUCCESS or SPASM_NO_SOLUTION
 */
int spasm_LU_solve(spasm * A, const spasm_ZZp * b, spasm_ZZp * x) {
	spasm_ZZp *y, *z, *w;
	spasm_lu *LU;
	spasm *L, *U;
	int n, m, r, i, ok;
	int *q, *p, *qinv;

	/* check inputs */
	assert(A != NULL);
	assert(b != NULL);
	n = A->n;
	m = A->m;

	qinv = spasm_malloc(m * sizeof(int));
	p = spasm_malloc(n * sizeof(int));
	spasm_find_pivots(A, p, qinv);
	LU = spasm_LU(A, p, SPASM_KEEP_L);
	L = LU->L;
	U = LU->U;

	r = U->n;

	/* get workspace */
	y = spasm_malloc(m * sizeof(spasm_ZZp));
	z = spasm_malloc(r * sizeof(spasm_ZZp));
	w = spasm_malloc(n * sizeof(spasm_ZZp));
	q = spasm_malloc(m * sizeof(int));

	for (i = 0; i < m; i++) {
		if (LU->qinv[i] != -1) {
			q[LU->qinv[i]] = i;
		}
	}

	/* z.U = b  (if possible) */
	for (i = 0; i < m; i++) {
		y[i] = b[i];
	}
	ok = spasm_dense_forward_solve(U, y, z, q);

	if (ok == SPASM_SUCCESS) {
		/* y.LU = b */
		spasm_dense_back_solve(L, z, w, LU->p);
		spasm_ipvec(p, w, x, n);
	}
	free(y);
	free(z);
	free(w);
	free(q);
	free(p);
	free(qinv);
	spasm_free_LU(LU);
	return ok;
}
