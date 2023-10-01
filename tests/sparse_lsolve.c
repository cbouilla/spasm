#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"
#include "test_tools.h"
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
        // load lower-triangular matrix matrix
        struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
        struct spasm_csr *L = spasm_compress(T);
        spasm_triplet_free(T);
        int n = L->n;
        int m = L->m;

        assert(n >= m); // lower-trapezoidal
        assert(spasm_is_lower_triangular(L));

        // load RHS
        T = spasm_triplet_alloc(1, m, 10, prime, true);
        spasm_add_entry(T, 0, 0, 1);
        spasm_add_entry(T, 0, m / 2, 2);
        spasm_add_entry(T, 0, m - 1, 3);
        struct spasm_csr *B = spasm_compress(T);
        spasm_triplet_free(T);

        int *xi = spasm_malloc(3*m * sizeof(*xi));
        for (int j = 0; j < 3 * m; j++)
                xi[j] = 0;

        spasm_ZZp *x = malloc(n * sizeof(*x));
        spasm_ZZp *y = malloc(m * sizeof(*y));
        for (int j = 0; j < n; j++)
                x[j] = 0;
        for (int j = 0; j < m; j++)
                y[j] = 0;
        int *qinv = malloc(m * sizeof(int));
        for (int j = 0; j < m; j++)
                qinv[j] = j;
        spasm_sparse_triangular_solve(L, B, 0, xi, x, qinv);

        /* check solution */
        spasm_xApy(x, L, y);
        spasm_scatter(B, 0, -1, y);
        for (int i = 0; i < m; i++)
                if (y[i] != 0) {
                        printf("not ok - sparse triangular L-solve (y[%d] == %d)\n", i, y[i]);
                        exit(1);
                }
        printf("ok - sparse triangular L-solve\n");
        spasm_csr_free(L);
        spasm_csr_free(B);
        free(xi);
        free(x);
        free(y);
        exit(EXIT_SUCCESS);
}
