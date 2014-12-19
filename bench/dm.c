#include <assert.h>
#include <stdio.h>
#include "spasm.h"

spasm_partition *process_square_part(const spasm * B, int rx, int ry, int cx, int cy, int *p, int *q) {
    spasm *C;
    spasm_partition *SCC;
    int i, k, l;

    C = spasm_submatrix(B, rx, ry, cx, cy, SPASM_IGNORE_VALUES);
    SCC = spasm_strongly_connected_components(C);
    k = SCC->nr;

    for (i = 0; i < k; i++) {
        l = SCC->rr[i + 1] - SCC->rr[i];
        printf("          --> SCC of size %d x %d\n", l, l);
    }

    /* update permutations of B */
    spasm_range_pvec(p, rx, ry, SCC->p);
    spasm_range_pvec(q, cx, cy, SCC->q);

    spasm_csr_free(C);
    return SCC;
}

spasm_cc *process_rectangular_part(const spasm * B, int ra, int rb, int ca, int cb, int *p, int *q, const int *jmatch) {
    spasm *M, *MM;
    int n, m, CC_k, i, k, rx, ry, cx, cy, C_n, C_m;
    int *M_jmatch, *MM_jmatch;
    int *CC_qinv;
    spasm_partition *CC;
    spasm_cc *result;


    M = spasm_submatrix(B, ra, rb, ca, cb, SPASM_IGNORE_VALUES);
    M_jmatch = spasm_submatching(jmatch, ra, rb, ca, cb);

    n = M->n;
    m = M->m;

    printf("*) M (%d x %d -- %d nnz) : \n", n, m, spasm_nnz(M));

    /* --- connected components of M */
    CC = spasm_connected_components(M, NULL, M_jmatch, NULL);
    CC_k = CC->nr;

    printf("   *) %d connected components\n", CC_k);

    /* permute M to expose the connected components */
    CC_qinv = spasm_pinv(CC->q, m);
    MM = spasm_permute(M, CC->p, CC_qinv, SPASM_IGNORE_VALUES);
    free(CC_qinv);

    result = spasm_malloc(sizeof(spasm_cc));
    result->CC = CC;
    result->SCC = spasm_malloc(CC_k * sizeof(spasm_partition *));

    for (i = 0; i < CC_k; i++) {

        /* process i-th connected component of M */
        result->SCC[i] = NULL;
        C_n = CC->rr[i + 1] - CC->rr[i];
        C_m = CC->cc[i + 1] - CC->cc[i];

        if (C_n == 0 || C_m == 0) {
            continue;
        }
        /* extract the (square) perfectly-matched part */
        k = spasm_min(C_n, C_m);
        cx = CC->cc[i];
        ry = CC->rr[i + 1];
        if (C_n == C_m) {
            /* square case: the matching is perfect */
            rx = CC->rr[i];
            cy = CC->cc[i + 1];
        } else if (C_n < C_m) {
            /* horizontal case: matched columns are on the left */
            rx = CC->rr[i];
            cy = CC->cc[i] + k;
        } else {
            /* vertical case: matched rows are on the bottom */
            rx = CC->rr[i + 1] - k;
            cy = CC->cc[i + 1];
        }

        printf("      --> %d x %d\n", C_n, C_m);
        result->SCC[i] = process_square_part(MM, rx, ry, cx, cy, CC->p, CC->q);
    }

    /* update permutations of B */
    spasm_range_pvec(p, ra, rb, CC->p);
    spasm_range_pvec(q, ca, cb, CC->q);

    /* cleanup */
    free(M_jmatch);
    spasm_csr_free(MM);
    spasm_csr_free(M);

    return result;
}

spasm_dm *full_DM(const spasm * A) {
    spasm *B, *A_t;
    spasm_partition *DM;
    int n, m;
    int *p, *q, *rr, *cc;
    int *qinv;
    int *imatch, *jmatch, *Bjmatch;
    spasm_dm *result;

    assert(A != NULL);
    n = A->n;
    m = A->m;

    A_t = spasm_transpose(A, SPASM_IGNORE_VALUES);

    /* --- Maximum matching ------------------------------------------------- */
    jmatch = spasm_malloc(n * sizeof(int));
    imatch = spasm_malloc(m * sizeof(int));

    if (n < m) {
        spasm_maximum_matching(A, jmatch, imatch);
    } else {
        spasm_maximum_matching(A_t, imatch, jmatch);
    }

    /* --- coarse DM decomposition ----------------------------------------- */
    DM = spasm_dulmage_mendelson(A, A_t, jmatch, imatch);
    p = DM->p;
    q = DM->q;
    rr = DM->rr;
    cc = DM->cc;

    spasm_csr_free(A_t);
    free(imatch);

    printf("got coarse\n");

    result = spasm_malloc(sizeof(spasm_dm *));
    result->DM = DM;
    result->H = NULL;
    result->S = NULL;
    result->V = NULL;

    qinv = spasm_pinv(q, m);
    B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);

    printf("got B\n");

    /* optimization : this could be done in-place */
    Bjmatch = spasm_permute_row_matching(n, jmatch, p, qinv);
    printf("got Bjmatch\n");

    free(qinv);
    printf("freed qinv\n");

    free(jmatch);
    printf("freed jmatch\n");

    /* ------------------- H --------------------- */
    if (cc[2] != cc[0]) {
        result->H = process_rectangular_part(B, rr[0], rr[1], cc[0], cc[2], p, q, Bjmatch);
    }
    /* --------------- S ----------------------- */
    if (rr[2] != rr[1]) {
        result->S = process_rectangular_part(B, rr[1], rr[2], cc[2], cc[3], p, q, Bjmatch);
    }
    /* ------------------- V --------------------- */
    if (rr[4] != rr[2]) {
        result->V = process_rectangular_part(B, rr[2], rr[4], cc[3], cc[4], p, q, Bjmatch);
    }
    /* cleanup */
    spasm_csr_free(B);
    free(Bjmatch);

    return result;
}


int main() {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n, m, *qinv;

    T = spasm_load_sms(stdin, -1);
    A = spasm_compress(T);
    spasm_triplet_free(T);

    n = A->n;
    m = A->m;

    printf("A : %d x %d with %d nnz\n", n, m, spasm_nnz(A));

    x = full_DM(A);

    qinv = spasm_pinv(x->DM->q, m);
    B = spasm_permute(A, x->DM->p, qinv, SPASM_IGNORE_VALUES);
    free(qinv);

    /*
     * todo:recalculer le matching permut Ã ©(pour l 'affichage en couleur,
     * le matching suffit)
     */
    FILE *f = fopen("plop.ppm", "w");
    spasm_save_ppm(f, m, n, B, x->DM);
    fclose(f);

    spasm_csr_free(B);
    spasm_csr_free(A);
    return 0;
}
