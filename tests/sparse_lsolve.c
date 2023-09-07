#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

int main(int argc, char **argv)
{
        int prime = 32003;
        // load lower-triangular matrix matrix
        spasm_triplet *T = spasm_load_sms(stdin, prime);
        spasm *L = spasm_compress(T);
        spasm_triplet_free(T);
        int n = L->n;
        int m = L->m;

        assert(n >= m); // lower-trapezoidal
        assert(spasm_is_lower_triangular(L));

        // load RHS
        T = spasm_triplet_alloc(1, m, 10, 32003, prime);
        spasm_add_entry(T, 0, 0, 1);
        spasm_add_entry(T, 0, m / 2, 2);
        spasm_add_entry(T, 0, m - 1, 3);
        spasm *B = spasm_compress(T);
        spasm_triplet_free(T);

        int *xi = spasm_malloc(3*m * sizeof(*xi));
        spasm_vector_zero(xi, 3*m);

        spasm_GFp *x = malloc(n * sizeof(*x));
        spasm_GFp *y = malloc(m * sizeof(*y));
        spasm_vector_zero(x, n);
        spasm_vector_zero(y, m);
        int *qinv = malloc(m * sizeof(int));
        for (int j = 0; j < m; j++)
                qinv[j] = j;
        spasm_sparse_triangular_solve(L, B, 0, xi, x, qinv);

        /* check solution */
        spasm_xApy(x, L, y);
        spasm_scatter(B->j, B->x, B->p[0], B->p[1], prime - 1, y, prime);
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
