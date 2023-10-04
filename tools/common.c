#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <err.h>

#include "spasm.h"
#include "common.h"

const char *argp_program_version = SPASM_VERSION;
const char *argp_program_bug_address = SPASM_BUG_ADDRESS;

struct argp_option input_options[] = {
	{0,               0,   0,     0, "Input problem", 1 },
	{"matrix",       'm', "FILE", 0, "Read the input matrix from FILE", 1 },
	{"modulus",      'p', "P",    0, "Perform arithmetic modulo P", 1 },
	{0}
};

/* Parse a single option --- rank program. */
error_t parse_input_opt(int key, char *arg, struct argp_state *state)
{
	struct input_matrix *arguments = state->input;
	switch (key) {
	case ARGP_KEY_INIT:  /* set defaults */
		arguments->filename = NULL;
		arguments->prime = 42013;
		break;
	case 'm':
		arguments->filename = arg;
		break;
	case 'p':
		arguments->prime = atoll(arg);
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

struct argp input_argp = { input_options, parse_input_opt, NULL, NULL, NULL, NULL, NULL};


/* The options of the echelonization code */
enum ech_opt_key {
	NO_LOW_RANK, NO_DENSE, NO_GPLU,
	MAX_ITER, DENSE_THR, MIN_PIV_RATIO,
	DENSE_BLKSZ, MIN_RANK_RATIO, MAX_ASPECT_RATIO
};

struct argp_option echelonize_options[] = {
	{0,                     0,                 0, 0, "Echelonization sub-algorithms", -2},
	{"no-low-rank-mode",    NO_LOW_RANK,       0, 0, "Disable the (dense) low-rank mode", -2 },
	{"no-dense-mode",       NO_DENSE,          0, 0, "Don't use FFPACK", -2 },
	{"no-GPLU",             NO_GPLU,           0, 0, "Don't use GPLU", -2 },

	{0,                     0,                 0,  0, "Main echelonization options", -3 },
	{"max-iterations",      MAX_ITER,         "N", 0, "Compute at most N sparse Schur complements ", -3},
	{"dense-threshold",     DENSE_THR,        "D", 0, "Use FFPACK when the density is greater than D", -3},
	{"min-pivot-proportion", MIN_PIV_RATIO,   "P", 0, "Stop if the proportion of pivots found is less than P", -3},

	{0,                     0,                 0,  0, "Dense code options", -4},
	{"dense-block-size",    DENSE_BLKSZ,      "N", 0, "Use dense matrices of at most N rows", -4},
	{"min-rank-ratio",      MIN_RANK_RATIO,   "R", 0, "Use low-rank mode if k rows have rank <= k * X", -4},
	{"max-aspect-ratio",    MAX_ASPECT_RATIO, "R", 0, "Use low-rank mode if #rows / #columns >= R", -4},
	
	{ 0 }
};

/* Parse a single option --- rank program. */
error_t parse_ech_opt(int key, char *arg, struct argp_state *state)
{
	struct echelonize_opts *opts = state->input;
	switch (key) {
	case ARGP_KEY_INIT:  /* set defaults */
		spasm_echelonize_init_opts(opts);
		break;
	case NO_LOW_RANK:
		opts->enable_tall_and_skinny = 0;
		break;
	case NO_DENSE:
		opts->enable_dense = 0;
		break;
	case NO_GPLU:
		opts->enable_GPLU = 0;
		break;
	case MAX_ITER:
		opts->max_round = atoi(arg);
		break;
	case DENSE_THR:
		opts->sparsity_threshold = atof(arg);
		break;
	case MIN_PIV_RATIO:
		opts->min_pivot_proportion = atof(arg);
		break;
	case DENSE_BLKSZ:
		opts->dense_block_size = atoi(arg);
		break;
	case MIN_RANK_RATIO:
		opts->low_rank_ratio = atof(arg);
		break;
	case MAX_ASPECT_RATIO:
		opts->tall_and_skinny_ratio = atof(arg);
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

struct argp echelonize_argp = { echelonize_options, parse_ech_opt, NULL, NULL, NULL, NULL, NULL};


struct spasm_triplet * load_input_matrix(struct input_matrix *in, u8 *hash)
{
	if (in->filename == NULL)
		return spasm_triplet_load(stdin, in->prime, hash);    /* load from stdin */
	FILE *f = fopen(in->filename, "r");
	if (f == NULL)
		err(1, "Cannot open %s", in->filename);
	struct spasm_triplet *T = spasm_triplet_load(f, in->prime, hash);
	fclose(f);
	return T;
}


FILE * open_input(const char *filename)
{
	if (filename == NULL)
		return stdin;
	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "Cannot open %s", filename);
	return f;
}


FILE * open_output(const char *filename)
{
	if (filename == NULL)
		return stdout;
	FILE *f = fopen(filename, "w");
	if (f == NULL)
		err(1, "Cannot open %s", filename);
	return f;
}