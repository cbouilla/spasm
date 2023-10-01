#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

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
        struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
        struct spasm_csr *A = spasm_compress(T);
        spasm_triplet_free(T);
  
        int n = A->n;
        int m = A->m;
        i64 *Ap = A->p;
        int *Aj = A->j;
        spasm_ZZp *Ax = A->x;
        spasm_datatype datatype = spasm_datatype_choose(prime);
        printf("# prime=%" PRId64 ", using datatype %s\n", prime, spasm_datatype_name(datatype));

        void *M = spasm_malloc(m * n * spasm_datatype_size(datatype));
        for (int i = 0; i < n * m; i++)
                spasm_datatype_write(M, i, datatype, 0);
        for (int i = 0; i < n*m; i++)
                assert(spasm_datatype_read(M, i, datatype) == 0);

        /* scatter A into M */
        for (int i = 0; i < n; i++)
                for (i64 k = Ap[i]; k < Ap[i + 1]; k++) {
                        int j = Aj[k];
                        spasm_datatype_write(M, i * m + j, datatype, Ax[k]);
                        // M[i * m + j] = Ax[k];
                }

        for (int i = 0; i < n; i++) {
                printf("# ");
                for (int j = 0; j < m; j++)
                        printf("%6d ", spasm_datatype_read(M, i * m + j, datatype));
                printf("\n");
        }

        size_t *p = spasm_malloc(m * sizeof(*p));
        size_t *qinv = spasm_malloc(m * sizeof(*qinv));
        int rank = spasm_ffpack_rref(prime, n, m, M, m, datatype, qinv);
        printf("# echelonized ; rank = %d\n", rank);

        /* dump output */
        for (int i = 0; i < n; i++) {
                printf("# ");
                for (int j = 0; j < m; j++)
                        printf("%6d ", spasm_datatype_read(M, i * m + j, datatype));
                printf("\n");
        }

        // for (int j = 0; j < n; j++)
        //         printf("# P[%d] = %d\n", j, P[j]);
        for (int j = 0; j < m; j++)
                printf("# Qt[%d] = %zd\n", j, qinv[j]);

        /* check that all rows of the input matrix belong to the row-space of U */  
        spasm_ZZp *x = spasm_malloc(m * sizeof(*x));
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
                                spasm_ZZp Mkl = spasm_datatype_read(M, k * m + l, datatype);
                                x[j] = spasm_ZZp_axpy(A->field, -alpha, Mkl, x[j]);
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
