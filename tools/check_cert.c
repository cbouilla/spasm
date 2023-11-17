#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <argp.h>

#include "spasm.h"
#include "common.h"

/* Program documentation. */
char doc[] = "Verify a rank certificate";

struct cmdline_args {
	/* input problem */
	struct input_matrix input;
	
	/* options specific to the rank program */
	char *cert_file;
	};

/* The options we understand. */
struct argp_option options[] = {
	{0,               0,   0,     0, "Rank options", 2 },
	{"certificate",  'c', "FILE", 0, "Read the rank certificate from FILE", 2 },
	{ 0 }
};

/* Parse a single option --- rank program. */
error_t parse_check_opt(int key, char *arg, struct argp_state *state)
{
	struct cmdline_args *arguments = state->input;
	switch (key) {
	case 'c':
		arguments->cert_file = arg;
		break;
	case ARGP_KEY_ARG:
		fprintf(stderr, "ERROR: invalid argument ``%s''\n", arg);
		exit(1);
	case ARGP_KEY_INIT:
		arguments->cert_file = NULL;
		state->child_inputs[0] = &arguments->input;
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

struct argp_child children_parsers[] = {
	{ &input_argp,      0, 0, 0 },
	{0}
};

/* argp parser for the rank program */
struct argp argp = { options, parse_check_opt, NULL, doc, children_parsers, NULL, NULL };


/** Computes the rank of the input matrix using the hybrid strategy described in [PASCO'17] */
int main(int argc, char **argv)
{
	/* process command-line options */
	struct cmdline_args args;
	argp_parse(&argp, argc, argv, 0, 0, &args);

	/* load input matrix */
	u8 hash[32];
	struct spasm_triplet *T = load_input_matrix(&args.input, hash);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	
	/* load certificate */
	FILE *f = open_input(args.cert_file);
	struct spasm_rank_certificate proof;
	spasm_rank_certificate_load(f, &proof);

	bool correct = spasm_certificate_rank_verify(A, hash, &proof);
	
	if (!correct)
		fprintf(stderr, "CORRECT certificate\n");
	else
		fprintf(stderr, "INCORRECT certificate\n");	
	return correct;
}