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
  
        /* build left kernel basis of A */
        spasm *At = spasm_transpose(A, SPASM_WITH_NUMERICAL_VALUES);
        int *qinv = spasm_malloc(n * sizeof(*qinv));
        struct echelonize_opts opts;
        spasm_echelonize_init_opts(&opts);
        opts.enable_tall_and_skinny = 0;
        opts.enable_dense = 0;
        spasm *U = spasm_echelonize(At, qinv, &opts);
        spasm *Ut = spasm_transpose(U, SPASM_WITH_NUMERICAL_VALUES);
        spasm *K = spasm_kernel(U, qinv);
        
        // S*At == M*U, donc A*St == Ut * Mt
        // Si x*Ut == 0, alors x*A*St == 0. Donc ce n'est pas normal que ça échoue.
        // Si x*A  == 0, alors x*Ut * Mt == 0, donc ou bien x*Ut == 0, ou bien (x*Ut) * Mt == 0
        // cette dernière chose est assez improbable

        /* rows of K form a basis of the right-kernel of At, hence the left kernel of A */
        
        /* test that they are really kernel vectors */
        spasm_GFp *x = spasm_malloc(n * sizeof(*x));
        spasm_GFp *y = spasm_malloc(m * sizeof(*y));
        const i64 *Kp = K->p;
        const int *Kj = K->j;
        const spasm_GFp *Kx = K->x;
        for (int i = 0; i < K->n; i++) {
                printf("# testing vector %d\n", i);

                /* scatter K[i] into x */
                spasm_vector_zero(x, n);
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

                // A  : n x m
                // At : m x n
                // U  : r x n      rowspan(U) == rowspan(At) ----> colspan(Ut) == colspan(A)
                // Ut : n x r      

                /* y <-- x.Ut */
                spasm_vector_zero(y, m);
                spasm_xApy(x, Ut, y);
                for (int i = 0; i < m; i++) {
                        if (y[i] != 0) {
                                printf("not ok - vector not in kernel (product with Ut, y[%d] = %d)\n", i, y[i]);
                                exit(1);
                        }
                }

                /* y <-- x.A */
                spasm_vector_zero(y, m);
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
