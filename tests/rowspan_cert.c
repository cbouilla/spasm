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
	spasm_triplet *T = spasm_load_sms(stdin, prime, hash);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	opts.L = 1;
	spasm_lu *fact = spasm_echelonize(A, &opts);
	assert(fact->L != NULL);
	
	spasm_rowspan_certificate *proof = spasm_certificate_rowspan_create(A, fact, 1337);
	bool correct = spasm_certificate_rowspan_verify(A, fact->U, proof);
	assert(correct);
}
