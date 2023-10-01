#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"
#include "common.h"

/* Program documentation. */
char doc[] = "Compute a kernel basis of a sparse matrix";

struct cmdline_args {
	/* input problem */
	struct input_matrix input;
	
	/* common options to all programs that call spasm_echelonize() */
	struct echelonize_opts opts;

	/* options specific to the kernel program */
	bool left;
	char *output_filename;
};

/* The options we understand. */
struct argp_option options[] = {
	{0,               0,  0,      0, "Kernel options", 2 },
	{"left",         'l', 0,      0, "Compute the left-kernel", 2},
	{"output",       'o', "FILE", 0, "Write the kernel basis in FILE", 2 },
	{ 0 }
};

/* Parse a single option --- rank program. */
error_t parse_ker_opt(int key, char *arg, struct argp_state *state)
{
	struct cmdline_args *arguments = state->input;
	switch (key) {
	case 'l':
		arguments->left = 1;
		break;
	case 'o':
		arguments->output_filename = arg;
		break;
	case ARGP_KEY_INIT:
		arguments->left = 0;
		arguments->output_filename = NULL;
		state->child_inputs[0] = &arguments->input;
		state->child_inputs[1] = &arguments->opts;
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

struct argp_child children_parsers[] = {
	{ &input_argp,      0, 0, 0 },
	{ &echelonize_argp, 0, 0, 0 },
	{ 0 }
};

/* argp parser for the rank program */
struct argp argp = { options, parse_ker_opt, NULL, doc, children_parsers, NULL, NULL };


int main(int argc, char **argv)
{
	/* process command-line options */
	struct cmdline_args args;
	argp_parse(&argp, argc, argv, 0, 0, &args);

	/* load input matrix */
	struct spasm_triplet *T = load_input_matrix(&args.input, NULL);
	if (args.left) {
		fprintf(stderr, "Left-kernel, transposing\n");
		spasm_triplet_transpose(T);
	}
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);

	/* echelonize A */
	struct spasm_lu *fact = spasm_echelonize(A, &args.opts);
	spasm_csr_free(A);

	/* kernel basis */
	struct spasm_csr *K = spasm_kernel(fact);
	fprintf(stderr, "Kernel basis matrix is %d x %d with %" PRId64 " nz\n", K->n, K->m, spasm_nnz(K));
	
	FILE *f = open_output(args.output_filename);
	spasm_csr_save(K, f);
}