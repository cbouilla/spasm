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
        spasm_dense_lu *LU = spasm_dense_LU_alloc(m, prime);
        spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));

        /* compute a dense "LU" factorisation of the input matrix */
        for (int i = 0; i < n; i++) {
                spasm_vector_zero(x, m);
                for (int px = Ap[i]; px < Ap[i+1]; px++)
                        x[Aj[px]] = Ax[px];
                spasm_dense_LU_process(LU, x);
        }

        /* check that all rows of the input matrix belong to the row-space of U */  
        for(int i = 0; i < n; i++) {
                spasm_vector_zero(x, m);
                for(int px = Ap[i]; px < Ap[i+1]; px++)
                        x[ Aj[px] ] = Ax[px];
                if (spasm_dense_LU_process(LU, x)) {
                        printf("not ok - rowspan(U) == rowspan(A) (row %d)\n", i);
                        exit(1);
                }
        }

        printf("ok - rowspan(U) == rowspan(A)\n");
        spasm_csr_free(A);
        spasm_dense_LU_free(LU);
        free(x);
        exit(EXIT_SUCCESS);
}
