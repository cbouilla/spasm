#include <assert.h>
#include <stdio.h>
#include "spasm.h"

void process_square_part(const spasm *B, int rx, int ry, int cx, int cy, int *p, int *q, int *jmatch) {
  spasm *C;
  spasm_partition *SCC;
  int i, k, l;

  C = spasm_submatrix(B, rx, ry, cx, cy, SPASM_IGNORE_VALUES);
  SCC = spasm_strongly_connected_components(C);
  k = SCC->nr;

    for(i = 0; i < k; i++) {
      l = SCC->rr[i + 1] - SCC->rr[i];
      printf("          --> SCC of size %d x %d\n", l, l);
    }

    /* update permutations of B */
    spasm_range_pvec(p, rx, ry, SCC->p);
    spasm_range_pvec(q, cx, cy, SCC->q);

    spasm_partition_free(SCC);
    spasm_csr_free(C);
}

void process_rectangular_part(const spasm *B, int ra, int rb, int ca, int cb, int *p, int *q, int *jmatch) {
  spasm *M, *MM;
  int n, m, CC_k, i, k, rx, ry, cx, cy, C_n, C_m;
  int *M_jmatch, *MM_jmatch;
  int  *CC_qinv;
  spasm_partition *CC;

  M = spasm_submatrix(B, ra, rb, ca, cb, SPASM_IGNORE_VALUES);
  M_jmatch = spasm_submatching(jmatch, ra, rb, ca, cb);

  n = M->n;
  m = M->m;

  printf("*) M (%d x %d -- %d nnz) : \n",  n, m, spasm_nnz(M));

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

    /* process i-th connected component of M */

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

    printf("      --> %d x %d\n", C_n, C_m);
    process_square_part(MM, rx, ry, cx, cy, CC->p, CC->q, M_jmatch);
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
  spasm *A, *B, *A_t;
  spasm_partition *DM;
  int n, m, h, s, v;
  int *rr, *cc;
  int *p, *pinv, *q, *qinv;
  int *imatch, *jmatch, *Bjmatch, *Bimatch;
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

   spasm_csr_free(A_t);

   h = (cc[2] != 0);
   s = (rr[2] != rr[1]);
   v = (rr[4] != rr[2]);

   pinv = spasm_pinv(p, n);
   qinv = spasm_pinv(q, m);
   B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);
   Bjmatch = spasm_permute_row_matching(n, jmatch, p, qinv);
   Bimatch = spasm_permute_column_matching(m, imatch, pinv, q);
   free(pinv);
   free(qinv);

  /* ------------------- H --------------------- */
  if (h) {
    process_rectangular_part(B, rr[0], rr[1], cc[0], cc[2], p, q, Bjmatch);
  }

  /* --------------- S ----------------------- */
  if (s) {
    printf("*) S (%d x %d) : \n",  rr[2] - rr[1], cc[3] - cc[2]);
    process_square_part(B, rr[1], rr[2], cc[2], cc[3], p, q, Bjmatch);
  }

  /* ------------------- V --------------------- */
  if (v) {
    process_rectangular_part(B, rr[2], rr[4], cc[3], cc[4], p, q, Bjmatch);
  }

  /* save result */
  spasm_csr_free(B);

  qinv = spasm_pinv(q, m);
  B = spasm_permute(A, p, qinv, SPASM_IGNORE_VALUES);
  free(qinv);

  // todo : recalculer le matching permuté (pour l'affichage en couleur, le matching suffit)
  FILE *f = fopen("plop.ppm", "w");
  spasm_save_ppm(f, m, n, B, DM);
  fclose(f);

  free(imatch);
  free(jmatch);
  spasm_partition_free(DM);
  spasm_csr_free(B);
  spasm_csr_free(A);
  return 0;
}
