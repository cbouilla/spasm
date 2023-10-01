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
        parse_command_line_options(argc, argv);

        struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
        struct spasm_csr *A = spasm_compress(T);
        spasm_triplet_free(T);

        int n = A->n;
        int m = A->m;
  
        /* build left kernel basis of A */
        struct spasm_csr *At = spasm_transpose(A, true);
        struct echelonize_opts opts;
        spasm_echelonize_init_opts(&opts);
        struct spasm_lu *fact = spasm_echelonize(At, &opts);
        struct spasm_csr *U = fact->U;
        struct spasm_csr *Ut = spasm_transpose(U, true);
        struct spasm_csr *K = spasm_kernel(fact);

        /* rows of K form a basis of the right-kernel of At, hence the left kernel of A */
        
        /* test that they are really kernel vectors */
        spasm_ZZp *x = spasm_malloc(n * sizeof(*x));
        spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
        const i64 *Kp = K->p;
        const int *Kj = K->j;
        const spasm_ZZp *Kx = K->x;
        for (int i = 0; i < K->n; i++) {
                printf("# testing vector %d\n", i);

                /* scatter K[i] into x */
                for (int j = 0; j < n; j++)
                        x[j] = 0;
                int nonzero = 0;
                for (i64 p = Kp[i]; p < Kp[i + 1]; p++) {
                        int j = Kj[p];
                        x[j] = Kx[p];
                        nonzero += (Kx[p] != 0);
                }
                if (nonzero == 0) {
                        printf("not ok - zero vector in kernel\n");
                        exit(1);
                }

                assert(Ut->n == n);
                assert(Ut->m <= m);

                /* y <-- x.Ut */
                for (int j = 0; j < m; j++)
                        y[j] = 0;
                spasm_xApy(x, Ut, y);
                for (int i = 0; i < m; i++) {
                        if (y[i] != 0) {
                                printf("not ok - vector not in kernel (product with Ut, y[%d] = %d)\n", i, y[i]);
                                exit(1);
                        }
                }

                /* y <-- x.A */
                for (int j = 0; j < m; j++)
                        y[j] = 0;
                spasm_xApy(x, A, y);
                for (int i = 0; i < m; i++) {
                        if (y[i] != 0) {
                                printf("not ok - vector not in kernel (product with A, y[%d] = %d)\n", i, y[i]);
                                exit(1);
                        }
                }
        }

        printf("ok - (right-)kernel basis\n");
        spasm_csr_free(A);
        spasm_csr_free(K);
        free(x);
        free(y);
        exit(EXIT_SUCCESS);
}
