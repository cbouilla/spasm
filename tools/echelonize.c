#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"

i64 prime = 42013;
bool rref = 0;

struct echelonize_opts opts;

void parse_command_line_options(int argc, char **argv)
{
	struct option longopts[] = {
		{"modulus", required_argument, NULL, 'p'},
		{"rref", no_argument, NULL, 'r'},
		{"no-greedy-pivot-search", no_argument, NULL, 'g'},
		{"no-low-rank-mode", no_argument, NULL, 'l'},
		{"dense-block-size", required_argument, NULL, 'd'},
		{"low-rank-start-weight", required_argument, NULL, 'w'},
		{"sparsity-threshold", required_argument, NULL, 's'},
		{NULL, 0, NULL, 0}
	};
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'p':
			prime = atoi(optarg);
			break;
		case 'r':
			rref = 1;
			break;
		case 'g':
			opts.enable_greedy_pivot_search = 0;
			break;
		case 'd':
			opts.dense_block_size = atoi(optarg);
			break;
		case 'l':
			opts.enable_tall_and_skinny = 0;
			break;
		case 's':
			opts.sparsity_threshold = atof(optarg);
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

	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	int m = A->m;

	/* echelonize A */
	struct spasm_lu *fact = spasm_echelonize(A, &opts);
	spasm_csr_free(A);

	if (rref) {
		/* compute the RREF */
		int *Rqinv = spasm_malloc(m * sizeof(int));
		struct spasm_csr *R = spasm_rref(fact, Rqinv);
		spasm_csr_save(R, stdout);
		spasm_csr_free(R);
		free(Rqinv);
	} else {
		spasm_csr_save(fact->U, stdout);
	}
	spasm_lu_free(fact);
	exit(EXIT_SUCCESS);
}
