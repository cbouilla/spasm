#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "spasm.h"

int main(int argc, char **argv) {
        spasm_triplet *T = spasm_load_sms(stdin, 42013);
        spasm *A = spasm_compress(T);
        spasm_triplet_free(T);
        int n = A->n;
        int m = A->m;
        i64 *Ap = A->p;
        int *Aj = A->j;
        int *p = spasm_malloc(n * sizeof(*p));
        int *qinv = spasm_malloc(m * sizeof(*qinv));
  
        int size = spasm_maximum_matching(A, p, qinv);

        /* check that size is correct */
        int k = 0;
        for (int i = 0; i < n; i++)
                if (p[i] >= 0)
                        k += 1;
        if (k != size) {
                printf("not ok - maximum matching - size of row matching is wrong (%d vs %d)\n", size, k);
                exit(1);
        }

        k = 0;
        for (int j = 0; j < m; j++)
                if (qinv[j] >= 0)
                        k += 1;
        if (k != size) {
                printf("not ok - maximum matching - size of column matching is wrong\n");
                exit(1);
        }

        /* check consistency */
        for (int j = 0; j < m; j++)
                if (qinv[j] >= 0)
                        if (p[qinv[j]] != j) {
                                printf("not ok - maximum matching - row and column matchings are inconsistent\n");
                                exit(1);
                        }
  
        /* check that the entries actually exist */
        for (int i = 0; i < n; i++) {
                if (p[i] < 0)
                        continue;
                int found = 0;
                for (i64 px = Ap[i]; px < Ap[i+1]; px++)
                        if (Aj[px] == p[i]) {
                                found = 1;
                                break;
                        }
                if (!found) {
                        printf("not ok - maximum matching - matching entries do not exist\n");
                        exit(1);
                }
        }
        free(p);
        free(qinv);
        spasm_csr_free(A);
        printf("ok - maximum matching\n");
        return 0;
}
