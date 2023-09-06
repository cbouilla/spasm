#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv)
{
        spasm_triplet *T = spasm_load_sms(stdin, 42013);
        spasm *A = spasm_compress(T);
        spasm_triplet_free(T);

        int n = A->n;
        int m = A->m;

        /* compute the RREF of A */
        int *qinv = spasm_malloc(m * sizeof(int));
        spasm *U = spasm_echelonize(A, qinv, NULL);   /* NULL = default options */
        spasm *Ut = spasm_transpose(U, SPASM_WITH_NUMERICAL_VALUES);
        
        /* build kernel basis */
        spasm *K = spasm_kernel(Ut, qinv);
        
        /* rows of K form a basis of the left-kernel of At */
        assert(K->m == m);
        assert(K->n == m - U->n);
        spasm_csr_free(U);
        spasm_csr_free(Ut);
        free(qinv);
        
        /* test that they are really kernel vectors */
        spasm *At = spasm_transpose(A, SPASM_WITH_NUMERICAL_VALUES);
        spasm_GFp *x = spasm_malloc(m * sizeof(*x));
        spasm_GFp *y = spasm_malloc(n * sizeof(*y));
        const i64 *Kp = K->p;
        const int *Kj = K->j;
        const spasm_GFp *Kx = K->x;
        for (int i = 0; i < K->n; i++) {
                printf("# testing vector %d\n", i);
                spasm_vector_zero(x, m);
                spasm_vector_zero(y, n);

                /* check that vector is not zero */
                if (spasm_row_weight(K, i) == 0) {
                        printf("not ok - empty vector in kernel\n");
                        exit(1);
                }

                /* scatter K[i] into x */
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

                /* y <-- x.A */
                spasm_gaxpy(At, x, y);

                for (int i = 0; i < n; i++) {
                        if (y[i] != 0) {
                                printf("not ok - vector not in kernel\n");
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
