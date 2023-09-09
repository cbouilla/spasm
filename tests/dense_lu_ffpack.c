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
        i64 *Ap = A->p;
        int *Aj = A->j;
        spasm_GFp *Ax = A->x;

        double *M = spasm_malloc(m * n * sizeof(*M));
        spasm_ffpack_setzero(prime, n, m, M, m);
        for (int i = 0; i < n*m; i++)
                assert(M[i] == 0);

        /* scatter A into M */
        for (int i = 0; i < n; i++)
                for (i64 k = Ap[i]; k < Ap[i + 1]; k++) {
                        int j = Aj[k];
                        M[i * m + j] = Ax[k];
                }

        for (int i = 0; i < n; i++) {
                printf("# ");
                for (int j = 0; j < m; j++)
                        printf("%6.0f ", M[i * m + j]);
                printf("\n");
        }

        size_t *p = spasm_malloc(m * sizeof(*p));
        size_t *qinv = spasm_malloc(m * sizeof(*qinv));
        int rank = spasm_ffpack_echelonize(prime, n, m, M, m, qinv);
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
        for (int j = 0; j < m; j++)
                printf("# Qt[%d] = %zd\n", j, qinv[j]);

        /* check that all rows of the input matrix belong to the row-space of U */  
        spasm_GFp *x = spasm_malloc(m * sizeof(*x));
        for (int i = 0; i < n; i++) {
                // scatter A[i] into x
                for (int j = 0; j < m; j++)
                        x[j] = 0;
                for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
                        int j = Aj[px];
                        x[j] = Ax[px];
                }
                for (int k = 0; k < rank; k++) {
                        int j = qinv[k];        /* column with the pivot */
                        int alpha = x[j];
                        x[j] = 0;
                        for (int l = rank; l < m; l++) {
                                int j = qinv[l];
                                x[j] = (int) (x[j] - alpha * M[k * m + l]) % prime;
                        }
                }
                printf("# row %2d --> (", i);
                for (int j = 0; j < m; j++)
                        printf("%8d", x[j]);
                printf(")\n");
                for (int j = 0; j < m; j++)
                         assert(x[j] == 0);
        }
        printf("ok - rowspan(A) contained in rowspan(U)\n");
        // spasm_csr_free(A);
        free(M);
        // free(x);
        exit(EXIT_SUCCESS);
}
