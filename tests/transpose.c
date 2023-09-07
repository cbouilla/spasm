#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"


int main(int argc, char **argv)
{
        spasm_triplet *T = spasm_load_sms(stdin, 42013);
        spasm *A = spasm_compress(T);
        spasm_triplet_free(T);

        int n = A->n;
        int m = A->m;

        spasm *B = spasm_transpose(A, SPASM_WITH_NUMERICAL_VALUES);
        spasm *C = spasm_transpose(B, SPASM_WITH_NUMERICAL_VALUES);

        if (spasm_is_lower_triangular(A)) {
                printf("# A is lower-triangular.\n");
                printf("Check that A.T is upper-triangular\n");
                if (!spasm_is_upper_triangular(B)) {
                        printf("not ok - A.T not upper-triangular \n");
                        exit(1);
                }
                printf("Check that A.T.T is lower-triangular\n");
                if (!spasm_is_lower_triangular(C)) {
                        printf("not ok - A.T.T not lower-triangular \n");
                        exit(1);
                }
        }

        if (spasm_is_upper_triangular(A)) {
                printf("# A is upper-triangular.\n");
                printf("Check that A.T is lower-triangular\n");
                if (!spasm_is_lower_triangular(B)) {
                        printf("not ok - A.T not lower-triangular \n");
                        exit(1);
                }
                printf("Check that A.T.T is upper-triangular\n");
                if (!spasm_is_upper_triangular(C)) {
                        printf("not ok - A.T.T not upper-triangular \n");
                        exit(1);
                }
        }

        printf("# check that A.T.T == A\n");
        spasm_GFp *x = spasm_malloc(n * sizeof(*x));
        spasm_GFp *y = spasm_malloc(m * sizeof(*y));
        spasm_GFp *z = spasm_malloc(m * sizeof(*z));

        for (int i = 0; i < n; i++) {
                spasm_vector_zero(x, n);
                spasm_vector_zero(y, m);
                spasm_vector_zero(z, m);
                x[i] = 1;
                spasm_xApy(x, A, y);
                spasm_xApy(x, C, z);

        for (int j = 0; j < m; j++)
                if (y[j] != z[j]) {
	               printf("not ok - (A.T).T != A \n");
	               exit(0);
                }
        }
        printf("ok - transpose \n");
        free(x);
        free(y);
        free(z);
        spasm_csr_free(A);
        spasm_csr_free(B);
        spasm_csr_free(C);
        exit(EXIT_SUCCESS);
}