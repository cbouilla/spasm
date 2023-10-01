#include <stdio.h>
#include <stdlib.h>

#include "spasm.h"

/* permute the columns to flip the matrix */

int main()
{
	struct spasm_triplet *A = spasm_triplet_load(stdin, 42013, NULL);

	int m = A->m;
	i64 nz = A->nz;
	int *Aj = A->j;

	for (i64 k = 0; k < nz; k++)
		Aj[k] = m - Aj[k] - 1;

	spasm_triplet_save(A, stdout);
	return EXIT_SUCCESS;
}
