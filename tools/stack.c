#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "spasm.h"

/** permute a matrix so that all empty rows and columns are removed. */

int main() {
	struct spasm_triplet *A;
	int n, m, nz, u, v;
	int *p, *q, *Ai, *Aj;

	A = spasm_triplet_load(stdin, 42013, NULL);

	n = A->n;
	m = A->m;
	nz = A->nz;
	Ai = A->i;
	Aj = A->j;

	fprintf(stderr, "matrix has advertized dimension %d x %d... ", n, m);
	fflush(stderr);


	/* allocate result */
	p = spasm_malloc(n * sizeof(int));
	q = spasm_malloc(m * sizeof(int));
	for (int i = 0; i < n; i++)
		p[i] = 0;
	for (int j = 0; j < m; j++)
		q[j] = 0;
	
	/* mark non-empty rows/columns */
	for (int k = 0; k < nz; k++) {
		p[Ai[k]] = 1;
		q[Aj[k]] = 1;
	}

	/* prefix-sum: p[i] = sum(p[k], k in range(i)) */
	u = 0;
	for (int i = 0; i < n; i++)
		if (p[i] > 0)
			p[i] = u++;

	v = 0;
	for (int j = 0; j < m; j++)
		if (q[j] > 0)
			q[j] = v++;


	/* modify matrix */
	A->n = u;
	A->m = v;
	fprintf(stderr, "but is in fact %d x %d\n", u, v);

	for (int k = 0; k < nz; k++) {
		Ai[k] = p[Ai[k]];
		Aj[k] = q[Aj[k]];
	}

	spasm_triplet_save(A, stdout);
	spasm_triplet_free(A);
	free(p);
	free(q);
	return 0;
}
