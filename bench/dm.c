#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include "spasm.h"

#if 0
/** computes a Dulmage-Mendelson decomposition */

void show(const spasm * M, spasm_cc * Y) {
    int i, j, a, b, c, d, e, f, g, h, r;

    for (i = 0; i < Y->CC->nr; i++) {
        a = Y->CC->rr[i];
        b = Y->CC->cc[i];
        c = Y->CC->rr[i + 1];
        d = Y->CC->cc[i + 1];

        printf("   *) Connected component (%d x %d) --- (%d, %d) to (%d, %d)\n", c - a, d - b, a, b, c, d);
        if (Y->SCC[i] != NULL) {
            for (j = 0; j < Y->SCC[i]->nr; j++) {
                e = Y->SCC[i]->rr[j];
                f = Y->SCC[i]->cc[j];
                g = Y->SCC[i]->rr[j + 1];
                h = Y->SCC[i]->cc[j + 1];
                r = 0; // subrank(M, e, f, g, h);
                if (g - e > 1) {
                    nontrivial_diag_size += g - e;
                    nontrivial_diag_rank += r;
                }
                if (g - e > 1) {
                    trivial_diag_rank += 1;
                }
                printf("       *) SCC (%d x %d, deffect %d) --- (%d, %d) to (%d, %d)\n", g - e, h - f, spasm_min(g - e, h - f) - r, e, f, g, h);
            }
        }
    }
}
#endif

int main(int argc, char **argv) {
    int ch;

    /* options descriptor */
    struct option longopts[5] = {
        {"permuted", no_argument, NULL, 'p'},
        {"verbose", no_argument, NULL, 'v'},
        {"tabulated", no_argument, NULL, 't'},
        {"image", required_argument, NULL, 'i'},
        {NULL, 0, NULL, 0}
    };

  
    char mode = 'p';
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {
        case 'p':
        case 'v':
        case 't':
        case 'i':
            mode = ch;
            break;
        default:
            printf("Unknown option\n");
            exit(1);
        }
    }

    spasm_triplet *T = spasm_load_sms(stdin, 42013);
    spasm * A = spasm_compress(T);
    spasm_triplet_free(T);

    int n = A->n;
    int m = A->m;

    spasm_dm *DM = spasm_dulmage_mendelsohn(A);
    int *rr = DM->rr;
    int *cc = DM->cc;

    int * qinv = spasm_pinv(DM->q, m);
    spasm * B = spasm_permute(A, DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
    free(qinv);
    spasm_csr_free(A);

    /* terse output with just the size of the largest diagonal block */
    if (mode == 't') {
//        printf("%5d \t %5d \t %6d \t %6d \t %.1f \t %6d \t %.1f\n", n, m, spasm_nnz(A), i, 100.0 * i / spasm_min(n, m), j, 1.0 * i / j);
    }
    if (mode == 'v') {
        printf("structural rank = %d\n", rr[2] + cc[4] - cc[3]);
        int h_n = rr[1] - rr[0];
        int h_m = cc[2] - cc[0];
        if (h_n > 0 && h_m > 0)
            printf("*) H (%d x %d)\n", h_n, h_m);
        
        int s_n = rr[2] - rr[1];
        int s_m = cc[3] - cc[2];
        if (s_n > 0 && s_m > 0)
            printf("*) S (%d x %d) : \n", s_n, s_m);
            /* Do something with S */
        
        int v_n = rr[4] - rr[2];
        int v_m = cc[4] - cc[3];
        if (v_n > 0 && v_m > 0)
            printf("*) V (%d x %d)\n", v_n, v_m);

    }
    if (mode == 'p') {
        spasm_save_csr(stdout, B);
    }
    /*if (img_file != NULL) {
        FILE *f = fopen(img_file, "w");
        spasm_save_ppm(f, B, x);
        fclose(f);
    }*/
    spasm_csr_free(B);

    return 0;
}
