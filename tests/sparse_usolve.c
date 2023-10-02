#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"
#include "test_tools.h"

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
        // load upper-triangular matrix matrix
        struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
        struct spasm_csr *U = spasm_compress(T);
        spasm_triplet_free(T);
        int n = U->n;
        int m = U->m;
        if (m == 0) {
                printf("SKIP --- empty matrix / useless\n");
                exit(EXIT_SUCCESS);
        }

        assert(n <= m); // upper-trapezoidal
        assert(spasm_is_upper_triangular(U));

        // load RHS
        T = spasm_triplet_alloc(1, m, 10, prime, true);
        spasm_add_entry(T, 0, 0, 1);
        spasm_add_entry(T, 0, m / 2, 2);
        spasm_add_entry(T, 0, m - 1, 3);
        struct spasm_csr *B = spasm_compress(T);
        spasm_triplet_free(T);

        int *xi = spasm_malloc(3*m * sizeof(*xi));
        for (int j = 0; j < 3*m; j++)
                xi[j] = 0;
        spasm_ZZp *x = malloc(m * sizeof(*x));
        spasm_ZZp *y = malloc(m * sizeof(*y));
        for (int j = 0; j < m; j++) {
                x[j] = 0;
                y[j] = 0;
        }

        int *qinv = malloc(m * sizeof(int));
        for (int j = 0; j < n; j++)
                qinv[j] = j;
        for (int j = n; j < m; j++)
                qinv[j] = -1;
        spasm_sparse_triangular_solve(U, B, 0, xi, x, qinv);

        /* check solution */
        spasm_xApy(x, U, y);
        for (int j = n; j < m; j++)
                y[j] = spasm_ZZp_add(B->field, y[j], x[j]);

        spasm_scatter(B, 0, -1, y);

        for (int i = 0; i < m; i++)
                if (y[i] != 0) {
                        printf("not ok - sparse triangular U-solve\n");
                        exit(1);
                }
        printf("ok - sparse triangular U-solve\n");
        spasm_csr_free(U);
        spasm_csr_free(B);
        free(xi);
        free(x);
        free(y);
        exit(EXIT_SUCCESS);
}
