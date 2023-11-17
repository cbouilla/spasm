#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <argp.h>

#include "spasm.h"
#include "common.h"

/* Program documentation. */
char doc[] = "Compute the rank of a sparse matrix";

struct cmdline_args {
	/* input problem */
	struct input_matrix input;
	
	/* common options to all programs that call spasm_echelonize() */
	struct echelonize_opts opts;

	/* options specific to the rank program */
	bool certificate;
	char *cert_file;           /* for the certificate */
	bool allow_transpose;	   /* transpose ON by default */
};

/* The options we understand. */
struct argp_option options[] = {
	{0,               0,   0,     0, "Rank options", 2 },
	{"no-transpose", 't', 0,      0, "Do not transpose the input matrix", 2},
	{"certificate",  'c', 0,      0, "Output a rank certificate", 2 },
	{"output",       'o', "FILE", 0, "Write the rank certificate in FILE", 2 },
	{ 0 }
};

/* Parse a single option --- rank program. */
error_t parse_rank_opt(int key, char *arg, struct argp_state *state)
{
	struct cmdline_args *arguments = state->input;
	switch (key) {
	case 't':
		arguments->allow_transpose = 0;
		break;
	case 'c':
		arguments->certificate = 1;
		break;
	case 'o':
		arguments->cert_file = arg;
		break;
	case ARGP_KEY_ARG:
		fprintf(stderr, "ERROR: invalid argument ``%s''\n", arg);
		exit(1);
	case ARGP_KEY_INIT:
		arguments->allow_transpose = 1;
		arguments->certificate = 0;
		arguments->cert_file = NULL;
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
struct argp argp = { options, parse_rank_opt, NULL, doc, children_parsers, NULL, NULL };


/** Computes the rank of the input matrix using the hybrid strategy described in [PASCO'17] */
int main(int argc, char **argv)
{
	/* process command-line options */
	struct cmdline_args args;
	argp_parse(&argp, argc, argv, 0, 0, &args);

	/* load input matrix */
	u8 hash[32];
	struct spasm_triplet *T = load_input_matrix(&args.input, hash);
	if (args.allow_transpose && (T->n < T->m)) {
		fprintf(stderr, "[rank] transposing matrix : ");
		fflush(stderr);
		spasm_triplet_transpose(T);
	}
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int m = A->m;
	char hnnz[8];
	spasm_human_format(spasm_nnz(A), hnnz);
	fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, hnnz);

	if (args.certificate)
		args.opts.L = 1;

	double start_time = spasm_wtime();
	struct spasm_lu *fact = spasm_echelonize(A, &args.opts);   /* NULL = default options */
	double end_time = spasm_wtime();
	fprintf(stderr, "done in %.3f s rank = %d\n", end_time - start_time, fact->U->n);
	

	if (args.certificate) {
		assert(spasm_factorization_verify(A, fact, 42));
		assert(spasm_factorization_verify(A, fact, 1337));
		assert(spasm_factorization_verify(A, fact, 21011984));
		fprintf(stderr, "generating certificate\n");
		struct spasm_rank_certificate *proof = spasm_certificate_rank_create(A, hash, fact);
		fprintf(stderr, "checking certificate\n");
		bool correct = spasm_certificate_rank_verify(A, hash, proof);
		if (correct)
			fprintf(stderr, "CORRECT certificate\n");
		else
			fprintf(stderr, "INCORRECT certificate\n");
		if (args.cert_file) {
			fprintf(stderr, "Saving certificate to %s\n", args.cert_file);
			FILE *out = open_output(args.cert_file);
			spasm_rank_certificate_save(proof, out);
		}
	}
	spasm_lu_free(fact);
	spasm_csr_free(A);
	return 0;
}
