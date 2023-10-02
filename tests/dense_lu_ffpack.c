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

        size_t *P = spasm_malloc(n * sizeof(*P));
        size_t *Qinv = spasm_malloc(m * sizeof(*Qinv));
        int r = spasm_ffpack_LU(prime, n, m, M, m, datatype, P, Qinv);
        printf("# echelonized ; rank = %d\n", r);

        /* dump output */
        for (int i = 0; i < n; i++) {
                printf("# ");
                for (int j = 0; j < m; j++)
                        printf("%6d ", spasm_datatype_read(M, i * m + j, datatype));
                printf("\n");
        }

        for (int j = 0; j < n; j++)
                printf("# P[%d] = %zd\n", j, P[j]);
        for (int j = 0; j < m; j++)
                printf("# Qt[%d] = %zd\n", j, Qinv[j]);

        /* build our LU factorization */
        struct spasm_csr *U = spasm_csr_alloc(r, n, n*m, prime, true);
        struct spasm_triplet *L = spasm_triplet_alloc(n, r, n*m, prime, true);
        i64 *Up = U->p;
        int *Li = L->i;
        int *Uj = U->j;
        int *Lj = L->j;
        i64 unz = 0;
        i64 lnz = 0;
        spasm_ZZp *Ux = U->x;
        spasm_ZZp *Lx = L->x;
        
        for (int i = 0; i < n; i++) {
                int pi = P[i];
                // int i = Pt[pi];
                // assert(Pt[inew] == i);
                for (int j = 0; j < spasm_min(i + 1, r); j++) {
                        spasm_ZZp Mij = spasm_datatype_read(M, i  * m + j, datatype);
                        if (Mij == 0)
                                continue;
                        assert(j < r || Mij == 0);
                        Li[lnz] = pi;
                        Lj[lnz] = j;
                        Lx[lnz] = Mij;
                        lnz += 1;
                }
        }
        L->nz = lnz;

        /* fill U */
        for (int i = 0; i < r; i++) {
                /* implicit 1 in U */
                Uj[unz] = Qinv[i];
                Ux[unz] = 1;
                unz += 1;
                for (int j = i+1; j < m; j++) {
                        int jnew = Qinv[j];
                        spasm_ZZp x = spasm_datatype_read(M, i * m + j, datatype);
                        Uj[unz] = jnew;
                        Ux[unz] = x;
                        unz += 1;
                }
                Up[i + 1] = unz;
        }
        struct spasm_lu fact;
        fact.U = U;
        fact.L = spasm_compress(L);
        fact.qinv = spasm_malloc(m * sizeof(*fact.qinv));
        fact.p = spasm_malloc(n * sizeof(*fact.p));

        for (int j = 0; j < n; j++)
                fact.p[j] = P[j];
        for (int j = 0; j < m; j++)
                fact.qinv[j] = Qinv[j];

        // assert(spasm_factorization_verify(A, &fact, 65537));
        spasm_ZZp *x = malloc(n * sizeof(*x));
        spasm_ZZp *y = malloc(m * sizeof(*y));
        spasm_ZZp *u = malloc(n * sizeof(*u));
        spasm_ZZp *v = malloc(m * sizeof(*v));

        /* check that A == L*U */
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                        x[j] = 0;
                        u[j] = 0;
                }
                for (int j = 0; j < m; j++) {
                        y[j] = 0;
                        v[j] = 0;
                }
                x[i] = 1;

                spasm_xApy(x, A, y);     // y <- x*A
                spasm_xApy(x, fact.L, u); // u <- x*L
                spasm_xApy(u, fact.U, v); // v <- (x*L)*U

                for (int j = 0; j < m; j++) {
                        // printf("# x*A[%4d] = %8d VS x*L[%4d] = %8d VS x*LU[%4d] = %8d\n", j, y[j], j, u[j], j, v[j]);
                        if (y[j] != v[j])
                                printf("\nmismatch on row %d, column %d \n", i, j);
                        assert(y[j] == v[j]);
                }
        }


        printf("ok - L*U == A\n");
        // spasm_csr_free(A);
        free(M);
        // free(x);
        exit(EXIT_SUCCESS);
}
