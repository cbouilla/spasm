#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv)
{
        spasm_triplet *T = spasm_load_sms(stdin, 42013);
        spasm *A = spasm_compress(T);
        spasm_triplet_free(T);
        int m = A->m;

        /* echelonize A */
        int *Uqinv = spasm_malloc(m * sizeof(int));
        spasm *U = spasm_echelonize(A, Uqinv, NULL);   /* NULL = default options */
        spasm_csr_free(A);

        /* compute the RREF */
        int *Rqinv = spasm_malloc(m * sizeof(int));
        spasm *R = spasm_rref(U, Uqinv, Rqinv);
        spasm_csr_free(U);
        free(Uqinv);
        
        /* build kernel basis from the RREF */
        spasm *K = spasm_kernel(R, Rqinv);
        char hnnz[8];
        spasm_human_format(spasm_nnz(K), hnnz);
        fprintf(stderr, "Saving kernel basis (%d x %d with %s nnz)\n", K->n, K->m, hnnz);
        spasm_save_csr(stdout, K);
        
        /* rows of K form a basis of the left-kernel of At */
        spasm_csr_free(K);
        spasm_csr_free(R);
        free(Rqinv);        
        exit(EXIT_SUCCESS);
}
