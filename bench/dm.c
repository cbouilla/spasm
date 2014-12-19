#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/* bad models :
332
175
205
457
*/

spasm_partition *process_square_part(const spasm * B, int rx, int ry, int cx, int cy, int *p, int *q) {
    spasm *C;
    spasm_partition *SCC;
    int i, k, l;

    C = spasm_submatrix(B, rx, ry, cx, cy, SPASM_IGNORE_VALUES);
    SCC = spasm_strongly_connected_components(C);
    k = SCC->nr;

    for (i = 0; i < k; i++) {
        l = SCC->rr[i + 1] - SCC->rr[i];
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

    /* --- connected components of M */
    CC = spasm_connected_components(M, NULL, M_jmatch, NULL);
    CC_k = CC->nr;

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

    result = spasm_malloc(sizeof(spasm_dm));
    result->DM = DM;
    result->H = NULL;
    result->S = NULL;
    result->V = NULL;

    qinv = spasm_pinv(q, m);
    B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);

    /* optimization : this could be done in-place */
    Bjmatch = spasm_permute_row_matching(n, jmatch, p, qinv);
    free(qinv);
    free(jmatch);

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


void show(spasm_cc *Y) {
  int i, j;

  for(i = 0; i < Y->CC->nr; i++) {
    int C_n = Y->CC->rr[i + 1] - Y->CC->rr[i];
    int C_m = Y->CC->cc[i + 1] - Y->CC->cc[i];
    printf("   *) Connected component (%d x %d)\n", C_n, C_m);
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	int SCC_n = Y->SCC[i]->rr[j + 1] - Y->SCC[i]->rr[j];
	int SCC_m = Y->SCC[i]->cc[j + 1] - Y->SCC[i]->cc[j];
	printf("       *) SCC (%d x %d)\n", SCC_n, SCC_m);
      }
    }
  }
}

int main() {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n, m, i, j, *qinv;

    T = spasm_load_sms(stdin, -1);
    A = spasm_compress(T);
    spasm_triplet_free(T);

    n = A->n;
    m = A->m;

    printf("A : %d x %d with %d nnz\n", n, m, spasm_nnz(A));

    x = full_DM(A);

    if (x->H != NULL) {
      int h_n = x->DM->rr[1] - x->DM->rr[0];
      int h_m = x->DM->cc[2] - x->DM->cc[0];
      printf("*) H (%d x %d) : \n", h_n, h_m);
      show(x->H);
    }
    if (x->S != NULL) {
      int s_n = x->DM->rr[2] - x->DM->rr[1];
      int s_m = x->DM->cc[3] - x->DM->cc[2];
      printf("*) S (%d x %d) : \n", s_n, s_m);
      show(x->S);
    }

    if (x->V != NULL) {
      int v_n = x->DM->rr[4] - x->DM->rr[2];
      int v_m = x->DM->cc[4] - x->DM->cc[3];
      printf("*) V (%d x %d) : \n", v_n, v_m);
      show(x->V);
    }


    qinv = spasm_pinv(x->DM->q, m);
    B = spasm_permute(A, x->DM->p, qinv, SPASM_IGNORE_VALUES);
    free(qinv);

    FILE *f = fopen("permuted.sms", "w");
    spasm_save_csr(f, B);
    fclose(f);

    /*
     * todo:recalculer le matching permut Ã ©(pour l 'affichage en couleur,
     * le matching suffit)
     */
    f = fopen("plop.ppm", "w");
    spasm_save_ppm(f, m, n, B, x->DM);
    fclose(f);

    spasm_csr_free(B);
    spasm_csr_free(A);
    return 0;
}
