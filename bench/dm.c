#include <assert.h>
#include <stdio.h>
#include "spasm.h"

void process_rectangular_part(const spasm *B, int ra, int rb, int ca, int cb, int *p, int *q, int *jmatch) {
  spasm *M, *MM, *C;
  int n, m, CC_k, SCC_k, i, j, k, l, rx, ry, cx, cy, C_n, C_m;
  int *M_jmatch, *MM_jmatch, *C_jmatch;
  int  *CC_qinv;
  spasm_partition *CC, *SCC;

  M = spasm_submatrix(B, ra, rb, ca, cb, SPASM_IGNORE_VALUES);
  M_jmatch = spasm_submatching(jmatch, ra, rb, ca, cb);

  n = M->n;
  m = M->m;

  printf("*) H (%d x %d -- %d nnz) : \n",  n, m, spasm_nnz(M));

  /* --- connected components of M */
  CC = spasm_connected_components(M, NULL, M_jmatch, NULL);
  CC_k = CC->nr;

  printf("   *) %d connected components\n", CC_k);

  /* permute M to expose the connected components */
  CC_qinv = spasm_pinv(CC->q, m);
  MM = spasm_permute(M, CC->p, CC_qinv, SPASM_IGNORE_VALUES);
  MM_jmatch = spasm_permute_row_matching(n, M_jmatch, CC->p, CC_qinv);
  free(CC_qinv);

  // TODO : update the matching

  for(i = 0; i < CC_k; i++) {

    /* process C_i, the i-th connected component of M */
    C_jmatch = spasm_submatching(MM_jmatch, CC->rr[i], CC->rr[i + 1], CC->cc[i], CC->cc[i + 1]);

    C_n = CC->rr[i + 1] - CC->rr[i];
    C_m = CC->cc[i + 1] - CC->cc[i];
    assert(C_n != C_m);

    /* extract the (square) perfectly-matched part */
    k = spasm_min(C_n, C_m);
    cx = CC->cc[i];
    ry = CC->rr[i + 1];
    if (C_n < C_m) {
      /* horizontal case: matched columns are on the left */
      rx = CC->rr[i];
      cy = CC->cc[i] + k;
    } else {
      /* vertical case: matched rows are on the bottom */
      rx = CC->rr[i + 1] - k;
      cy = CC->cc[i + 1];
    }
    C = spasm_submatrix(MM, rx, ry, cx, cy, SPASM_IGNORE_VALUES);

    SCC = spasm_strongly_connected_components(C);
    SCC_k = SCC->nr;
    printf("      --> %d x %d with %d SCC\n", C_n, C_m,  SCC_k);

    for(j = 0; j < SCC_k; j++) {
      l = SCC->rr[j + 1] - SCC->rr[j];
      printf("          --> %d x %d\n", l, l);
    }

    /* update permutations of M */
    spasm_range_pvec(CC->p, rx, ry, SCC->p);
    spasm_range_pvec(CC->q, cx, cy, SCC->q);

    spasm_partition_free(SCC);
    spasm_csr_free(C);
  }

  free(MM_jmatch);
  spasm_csr_free(MM);

  /* update permutations of B */
  spasm_range_pvec(p, ra, rb, CC->p);
  spasm_range_pvec(q, ca, cb, CC->q);


  /* TODO : update jmatch */
  free(M_jmatch);

  spasm_partition_free(CC);
  spasm_csr_free(M);
}

int main() {
  spasm_triplet *T;
  spasm *A, *B, *S, *V, *A_t;
  spasm_partition *DM, *P;
  int n, m, h, s, v, k, largest, ns, i, j;
  int *rr, *Prr, *cc, *Pcc;
  int *p, *pinv, *q, *qinv, *Vp, *Vq;
  int *imatch, *jmatch, *Bjmatch, *Bimatch,  *Vjmatch;
  double start;

  T = spasm_load_sms(stdin, -1);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  n = A->n;
  m = A->m;

  printf("A : %d x %d with %d nnz\n", n, m, spasm_nnz(A));

  A_t = spasm_transpose(A, SPASM_IGNORE_VALUES);

  /* --- Maximum matching ------------------------------------------------- */
  jmatch = spasm_malloc(n * sizeof(int));
  imatch = spasm_malloc(m * sizeof(int));
  start = spasm_wtime();

   if (n < m) {
     fprintf(stderr, "Maximum matching of A : ");
     fflush(stderr);
     spasm_maximum_matching(A, jmatch, imatch);
   } else {
     fprintf(stderr, "Maximum matching of A.T : ");
     fflush(stderr);
     spasm_maximum_matching(A_t, imatch, jmatch);
   }
   fprintf(stderr, "%.1f s\n", spasm_wtime() - start);

   /* --- coarse DM decomposition ----------------------------------------- */
   DM = spasm_dulmage_mendelson(A, A_t, jmatch, imatch);
   p = DM->p;
   q = DM->q;
   rr = DM->rr;
   cc = DM->cc;

   h = (cc[2] != 0);
   s = (rr[2] != rr[1]);
   v = (rr[4] != rr[2]);

   pinv = spasm_pinv(p, n);
   qinv = spasm_pinv(q, m);
   B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);
   Bjmatch = spasm_permute_row_matching(n, jmatch, p, qinv);
   Bimatch = spasm_permute_column_matching(m, imatch, pinv, q);


  /* ------------------- H --------------------- */
  if (h) {
    process_rectangular_part(B, rr[0], rr[1], cc[0], cc[2], p, q, Bjmatch);
  }

  /* --------------- S ----------------------- */
  if (s) {
    S = spasm_submatrix(B, rr[1], rr[2], cc[2], cc[3], SPASM_IGNORE_VALUES);
    printf("*) S (%d x %d --- %d nnz)\n", S->n, S->m, spasm_nnz(S));

    P = spasm_strongly_connected_components(S);
    Prr = P->rr;
    k = P->nr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, Prr[i + 1] - Prr[i]);
      ns += (Prr[i + 1] - Prr[i] > 1);
    }
    if (k > 1) {
      printf("   * %d strongly connected components, %d non-singleton, largest = %.1f %%\n", k, ns, 100.0 * largest / S->n);
      for(j = 0; j < k; j++) {
	printf("     --> %d x %d\n", Prr[j + 1] - Prr[j], Prr[j + 1] - Prr[j]);
      }

      spasm_range_pvec(p, rr[1], rr[2], P->p);
      spasm_range_pvec(q, cc[2], cc[3], P->q);
    }

    spasm_partition_free(P);
    spasm_csr_free(S);
  }

  /* ------------------- V --------------------- */
  if (v) {
    V = spasm_submatrix(B, rr[2], rr[4], cc[3], cc[4], SPASM_IGNORE_VALUES);
    printf("*) V (%d x %d -- %d nnz) : \n",  V->n, V->m, spasm_nnz(V));

    /* --- translate the maximum matching to a column-perfect matching of V */
    Vjmatch = spasm_submatching(Bjmatch, rr[2], rr[4], cc[3], cc[4]);

    P = spasm_connected_components(V, NULL, Vjmatch, NULL);
    Vp = P->p;
    Vq = P->q;
    Prr = P->rr;
    Pcc = P->cc;

    if (P->nr > 1) {
      printf("   *) %d connected components : \n", P->nr);
      for(i = 0; i < P->nr; i++) {
	printf("      --> %d x %d\n", Prr[i + 1] - Prr[i], Pcc[i + 1] - Pcc[i]);
      }
      /* in place permute p to reflect block decomposition */
      spasm_range_pvec(p, rr[2], rr[4], Vp);
      spasm_range_pvec(q, cc[3], rr[4], Vq);
    }
    spasm_csr_free(V);
    spasm_partition_free(P);
  }

  /* save result */
  spasm_csr_free(B);

  qinv = spasm_pinv(q, m);
  B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);
  spasm_csr_free(A);
  free(qinv);

  // todo : recalculer le matching permuté (pour l'affichage en couleur, le matching suffit)
  FILE *f = fopen("plop.ppm", "w");
  spasm_save_ppm(f, m, n, B, DM);
  fclose(f);

  spasm_partition_free(DM);
  spasm_csr_free(B);
  return 0;
}
