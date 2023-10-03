#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"

i64 prime = 42013;

void parse_command_line_options(int argc, char **argv)
{
        struct option longopts[] = {
                {"modulus", required_argument, NULL, 'p'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'p':
                        prime = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
}

int main(int argc, char **argv)
{
	u8 hash[32];
	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, hash);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);

	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	opts.L = 1;
	struct spasm_lu *fact = spasm_echelonize(A, &opts);
	assert(fact->L != NULL);

	struct spasm_rank_certificate *proof = spasm_certificate_rank_create(A, hash, fact);
	bool correct = spasm_certificate_rank_verify(A, hash, proof);
	assert(correct);
}
