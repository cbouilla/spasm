#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <argp.h>

#include "spasm.h"
#include "common.h"

/* Program documentation. */
char doc[] = "Solve sparse linear systems";

struct cmdline_args {
	/* input problem */
	struct input_matrix input;
	
	/* common options to all programs that call spasm_echelonize() */
	struct echelonize_opts opts;

	/* options specific to the solve program */
	char *rhs_filename;
	char *output_filename;
};

/* The options we understand. */
struct argp_option options[] = {
	{0,               0,   0,     0, "solve options", 2 },
	{"rhs",          'r', "FILE", 0, "Load the RHS matrix from FILE", 2 },
	{"output",       'o', "FILE", 0, "Write the solution matrix in FILE", 2 },
	{ 0 }
};

/* Parse a single option --- solve program. */
error_t parse_solve_opt(int key, char *arg, struct argp_state *state)
{
	struct cmdline_args *arguments = state->input;
	switch (key) {
	case 'r':
		arguments->rhs_filename = arg;
		break;
	case 'o':
		arguments->output_filename = arg;
		break;
	case ARGP_KEY_ARG:
		fprintf(stderr, "ERROR: invalid argument ``%s''\n", arg);
		exit(1);
	case ARGP_KEY_INIT:
		arguments->rhs_filename = NULL;
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
	{0}
};

/* argp parser for the rank program */
struct argp argp = { options, parse_solve_opt, NULL, doc, children_parsers, NULL, NULL };


/** solve several sparse linear systems */
int main(int argc, char **argv)
{
	/* process command-line options */
	struct cmdline_args args;
	argp_parse(&argp, argc, argv, 0, 0, &args);

	fprintf(stderr, "Loading A\n");
	u8 hash[32];
	struct spasm_triplet *T = load_input_matrix(&args.input, hash);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;

	fprintf(stderr, "Loading B\n");
	struct input_matrix rhs_in;
	rhs_in.prime = args.input.prime;
	rhs_in.filename = args.rhs_filename;
	T = load_input_matrix(&rhs_in, hash);
	struct spasm_csr *B = spasm_compress(T);
	spasm_triplet_free(T);

	/* echelonize A */
	fprintf(stderr, "Echelonizing A\n");
	char hnnz[8];
	spasm_human_format(spasm_nnz(A), hnnz);
	fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, hnnz);
	args.opts.L = 1;
	double start_time = spasm_wtime();
	struct spasm_lu *fact = spasm_echelonize(A, &args.opts);   /* NULL = default options */
	double end_time = spasm_wtime();
	fprintf(stderr, "echelonization done in %.3f s rank = %d\n", end_time - start_time, fact->U->n);
	
	fprintf(stderr, "Solving XA == B\n");
	bool *ok = spasm_malloc(B->n * sizeof(*ok));
	struct spasm_csr *X = spasm_gesv(fact, B, ok);
	for (int i = 0; i < B->n; i++)
		if (!ok[i])
			fprintf(stderr, "WARNING: no solution for row %d\n", i);
	fprintf(stderr, "done\n");

	FILE *f = open_output(args.output_filename);
	spasm_csr_save(X, f);

	exit(EXIT_SUCCESS);
}
