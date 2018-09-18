#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) 
{
	assert(argc > 1);
	int root = atoi(argv[1]);

	/*
	 * computes the list of reachable columns starting from column root. 
	 * pivots are assumed to be on the diagonal.
	 */

	spasm_triplet *T = spasm_load_sms(stdin, 42013);
	spasm *U = spasm_compress(T);
	spasm_triplet_free(T);

	int n = U->n;
	int m = U->m;
	assert(n <= m); /* why ? */

	int *pstack = malloc(m * sizeof(int));
	int *xj = malloc(m * sizeof(int));
	int *marks = malloc(m * sizeof(int));
	for (int j = 0 ; j < m; j++)
		marks[j] = 0;

	int *qinv = malloc(m * sizeof(int));
	for (int j = 0; j < n; j++)
		qinv[j] = j;
	for (int j = n; j < m; j++)
		qinv[j] = -1;

	int top = spasm_dfs(root, U, m, xj, pstack, marks, qinv);
	for( ; top < m; top++)
		printf("%d\n", xj[top]);

	spasm_csr_free(U);
	free(pstack);
	free(xj);
	free(marks);
	return 0;
}
