#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"

int prime = 42013;
bool rref = 0;

void parse_command_line_options(int argc, char **argv)
{
	struct option longopts[] = {
		{"modulus", required_argument, NULL, 'p'},
		{"rref", no_argument, NULL, 'r'},
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
		default:
			errx(1, "Unknown option\n");
		}
	}
}

int main(int argc, char **argv)
{
	parse_command_line_options(argc, argv);

        spasm_triplet *T = spasm_load_sms(stdin, prime);
        spasm *A = spasm_compress(T);
        spasm_triplet_free(T);
        int m = A->m;

        /* echelonize A */
        int *Uqinv = spasm_malloc(m * sizeof(int));
        spasm *U = spasm_echelonize(A, Uqinv, NULL);   /* NULL = default options */
        spasm_csr_free(A);

        if (rref) {
        	/* compute the RREF */
        	int *Rqinv = spasm_malloc(m * sizeof(int));
        	spasm *R = spasm_rref(U, Uqinv, Rqinv);
        	spasm_csr_free(U);
        	U = R;
        	free(Rqinv);        
        }
        free(Uqinv);
        exit(EXIT_SUCCESS);
}
