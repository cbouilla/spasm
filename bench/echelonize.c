#include <assert.h>
#include <stdio.h>
#include "spasm.h"
#include <getopt.h>

/** produces a matrix U such that:
    *) rowspan(U) == rowspan(A) 
    *) there exists a column permutation Q such that U.Q is upper-trapezoidal 
*/
int main(int argc, char **argv)
{
	char nnz[6];
	int n_times = -1;
	int prime = 42013;

	/* options descriptor */
	struct option longopts[3] = {
		{"max-recursion", required_argument, NULL, 'm'},
		{"modulus", required_argument, NULL, 'p'},
		{NULL, 0, NULL, 0}
	};

	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'm':
			n_times = atoi(optarg);
			break;
		case 'p':
			prime = atoi(optarg);
			break;
		default:
			printf("Unknown option\n");
			exit(1);
		}
	}

	spasm_triplet *T = spasm_load_sms(stdin, prime);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	int n = A->n;
	int m = A->m;
	spasm_human_format(spasm_nnz(A), nnz);
	fprintf(stderr, "A is %d x %d (%s nnz)\n", n, m, nnz);

	int *qinv = spasm_malloc(m * sizeof(*qinv));
	spasm_lu *LU = spasm_echelonize(A, n_times);
	free(qinv);

	spasm_save_csr(stdout, LU->U);
	spasm_free_LU(LU);
	return 0;
}