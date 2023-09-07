#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"

int prime = 42013;
struct echelonize_opts opts;
bool left = 0;

void parse_command_line_options(int argc, char **argv)
{
	struct option longopts[] = {
		{"modulus", required_argument, NULL, 'p'},
		{"left", no_argument, NULL, 'L'},
		{"no-greedy-pivot-search", no_argument, NULL, 'g'},
		{"no-low-rank-mode", no_argument, NULL, 'l'},
		{"dense-block-size", required_argument, NULL, 'd'},
		{"low-rank-start-weight", required_argument, NULL, 'w'},
		{NULL, 0, NULL, 0}
	};
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'p':
			prime = atoi(optarg);
			break;
		case 'L':
			left = 1;
			break;
		case 'g':
			opts.enable_greedy_pivot_search = 0;
			break;
		case 'd':
			opts.dense_block_size = atoi(optarg);
			fprintf(stderr, "Using dense block size %d\n", opts.dense_block_size);
			break;
		case 'l':
			opts.enable_tall_and_skinny = 0;
			break;
		case 'w':
			opts.low_rank_start_weight = atoi(optarg);
			break;
		default:
			errx(1, "Unknown option\n");
		}
	}
}

int main(int argc, char **argv)
{
	spasm_echelonize_init_opts(&opts);
	parse_command_line_options(argc, argv);

	spasm_triplet *T = spasm_load_sms(stdin, prime);
	if (left) {
		fprintf(stderr, "Left-kernel, transposing\n");
		spasm_triplet_transpose(T);
	}
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);
	int m = A->m;

	/* echelonize A */
	int *qinv = spasm_malloc(m * sizeof(*qinv));
	spasm *U = spasm_echelonize(A, qinv, &opts);
	spasm_csr_free(A);

        spasm *K = spasm_kernel(U, qinv);
	spasm_save_csr(stdout, K);
	fprintf(stderr, "Kernel basis matrix is %d x %d with %" PRId64 " nz\n", K->n, K->m, spasm_nnz(K));
	free(qinv);
	spasm_csr_free(U);
	spasm_csr_free(K);
        exit(EXIT_SUCCESS);
}
