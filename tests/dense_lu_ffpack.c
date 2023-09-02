#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv)
{
        int prime = 42013;
        spasm_triplet *T = spasm_load_sms(stdin, prime);
        spasm *A = spasm_compress(T);
        spasm_triplet_free(T);
  
        int n = A->n;
        int m = A->m;
        int *Ap = A->p;
        int *Aj = A->j;
        spasm_GFp *Ax = A->x;

        double *M = spasm_malloc(m * n * sizeof(*M));
        spasm_ffpack_setzero(prime, n, m, M, m);
        for (int i = 0; i < n*m; i++)
                assert(M[i] == 0);

        /* scatter A into M */
        for (int i = 0; i < n; i++)
                for (int k = Ap[i]; k < Ap[i + 1]; k++) {
                        int j = Aj[k];
                        M[i * m + j] = Ax[k];
                }

        // for (int i = 0; i < n; i++) {
        //         printf("# ");
        //         for (int j = 0; j < m; j++)
        //                 printf("%6.0f ", M[i * m + j]);
        //         printf("\n");
        // }

        size_t *Q = spasm_malloc(m * sizeof(*Q));
        int rank = spasm_ffpack_echelonize(prime, n, m, M, m, Q);
        printf("# echelonized ; rank = %d\n", rank);

        /* dump output */
        for (int i = 0; i < n; i++) {
                printf("# ");
                for (int j = 0; j < m; j++)
                        printf("%6.0f ", M[i * m + j]);
                printf("\n");
        }

        // for (int j = 0; j < n; j++)
        //         printf("# P[%d] = %d\n", j, P[j]);
        // for (int j = 0; j < m; j++)
        //         printf("# Q[%d] = %zd\n", j, Q[j]);

        /* check if really triangular */
        for (int i = 0; i < rank; i++) {
                int pivot_col = Q[i];
                // printf("# Q[%d] = %zd\n", i, Q[i]);
                assert(pivot_col >= i);
                // for (int j = 0; j < pivot_col; j++)
                //         assert(M[i*m + j] == 0);
                // assert(M[i*m + pivot_col] == 1);
        }

        /* check that all rows of the input matrix belong to the row-space of U */  
        spasm_GFp *x = spasm_malloc(m * sizeof(*x));
        for (int i = 0; i < n; i++) {
                spasm_vector_zero(x, m);
                for (int px = Ap[i]; px < Ap[i + 1]; px++) {
                        int j = Aj[px];
                        x[j] = Ax[px];
                }
                for (int k = 0; k < rank; k++) {
                        int j = Q[k];
                        int alpha = x[j];
                        x[j] = 0;
                        for (int l = j + 1; l < m; l++)
                                x[l] = (int) (x[l] - alpha * M[k * m + l]) % prime;
                }
                for (int j = 0; j < m; j++)
                        assert(x[j] == 0);
        }
        printf("ok - rowspan(A) contained in rowspan(U)\n");
        // spasm_csr_free(A);
        free(M);
        // free(x);
        exit(EXIT_SUCCESS);
}
