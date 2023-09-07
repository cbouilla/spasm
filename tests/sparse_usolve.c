#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

int main(int argc, char **argv)
{
        // load upper-triangular matrix matrix
        spasm_triplet *T = spasm_load_sms(stdin, 32003);
        spasm *U = spasm_compress(T);
        spasm_triplet_free(T);
        int n = U->n;
        int m = U->m;

        assert(n<= m); // upper-trapezoidal
        assert(spasm_is_upper_triangular(U));

        // load RHS
        T = spasm_triplet_alloc(1, m, 10, 32003, true);
        spasm_add_entry(T, 0, 0, 1);
        spasm_add_entry(T, 0, m / 2, 2);
        spasm_add_entry(T, 0, m - 1, 3);
        spasm *B = spasm_compress(T);
        spasm_triplet_free(T);

        int *xi = spasm_malloc(3*m * sizeof(*xi));
        spasm_vector_zero(xi, 3*m);

        spasm_GFp *x = malloc(m * sizeof(*x));
        spasm_GFp *y = malloc(m * sizeof(*y));
        spasm_vector_zero(x, m);
        spasm_vector_zero(y, m);
        int *qinv = malloc(m * sizeof(int));
        for (int j = 0; j < n; j++)
                qinv[j] = j;
        for (int j = n; j < m; j++)
                qinv[j] = -1;
        spasm_sparse_triangular_solve(U, B, 0, xi, x, qinv);

        /* check solution */
        spasm_xApy(x, U, y);
        for (int j = n; j < m; j++)
                y[j] = (y[j] + x[j]) % B->prime;

        spasm_scatter(B->j, B->x, B->p[0], B->p[1], B->prime - 1, y, B->prime);

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
