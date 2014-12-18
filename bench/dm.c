#include <assert.h>
#include <stdio.h>
#include "spasm.h"

int main() {
  spasm_triplet *T;
  spasm *A, *B, *H, *S, *V, *A_t, *HH, *C;
  spasm_partition *DM, *P;
  int n, m, h, s, v, k, largest, ns, i, j, l;
  int *rr, *Prr, *cc, *Pcc;
  int *p, *pinv, *q, *qinv, *Hp, *Hq, *Hqinv, *Vp, *Vq, *Sp, *Sq;
  int *imatch, *jmatch, *Bjmatch, *Bimatch, *Himatch, *Hjmatch, *HHjmatch, *Cjmatch;
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
    H = spasm_submatrix(B, rr[0], rr[1], cc[0], cc[2], SPASM_IGNORE_VALUES);
    printf("*) H (%d x %d -- %d nnz) : \n",  H->n, H->m, spasm_nnz(H));

    /* --- translate the maximum matching to a row-perfect matching of H */
    Hjmatch = spasm_submatching(Bjmatch, rr[0], rr[1]);

    /* --- connected components of H */
    P = spasm_connected_components(H, NULL, Hjmatch, NULL);
    Prr = P->rr;
    Pcc = P->cc;
    Hp = P->p;
    Hq = P->q;
    if (P->nr > 1) {
      printf("   *) %d connected components\n", P->nr);

      Hqinv = spasm_pinv(Hq, H->m);
      HH = spasm_permute(H, Hp, Hqinv, SPASM_IGNORE_VALUES);
      HHjmatch = spasm_permute_row_matching(H->n, Hjmatch, Hp, Hqinv);

      for(k = 0; k < P->nr; k++) {
	printf("      --> %d x %d\n", Prr[k + 1] - Prr[k], Pcc[k + 1] - Pcc[k]);

	Cjmatch = spasm_submatching(HHjmatch, Prr[k], Prr[k + 1]);
	printf("      --> column matched : ");
	for(i = Prr[k]; i < Prr[k + 1]; i++) {
	  printf("%d ", Cjmatch[ i - Prr[k] ]);
	}
	printf("\n");
	// size of the (square) column-matched part
	l = Prr[k + 1] - Prr[k];
	C = spasm_submatrix(HH, Prr[k], Prr[k+1], Pcc[k], Pcc[k] + l, SPASM_IGNORE_VALUES);
	printf("      --> size of perfectly matched part : %d x %d with %d NNZ\n", l, l, spasm_nnz(C));
	// calculer ses SCC
	// mettre à jour les permutations de H
	spasm_csr_free(C);
      }

      free(Hqinv);
      free(HHjmatch);
      spasm_csr_free(HH);
    }
    spasm_csr_free(H);
    spasm_partition_free(P);

    // une fois que tout est fini, mettre à jour p et q pour la matrice de départ A
    // FIXME
    /* in place permute p to reflect block decomposition */
      for(i = 0; i < H->n; i++) {
	Hp[i] = p[ Hp[i] ];
      }
      for(i = 0; i < H->n; i++) {
	p[i] = Hp[i];
      }

      /* in place permute q to reflect block decomposition */
      for(j = 0; j < H->m; j++) {
	Hq[j] = q[ Hq[j] ];
      }
      for(j = 0; j < H->m; j++) {
	q[j] = Hq[j];
      }
      
      free(Hjmatch);
  }

  /* --------------- S ----------------------- */
  if (s) {
    S = spasm_submatrix(B, rr[1], rr[2], cc[2], cc[3], SPASM_IGNORE_VALUES);
    printf("*) S (%d x %d --- %d nnz)\n", S->n, S->m, spasm_nnz(S));

    P = spasm_strongly_connected_components(S);
    Prr = P->rr;
    Sp = P->p;
    Sq = P->q;
    k = P->nr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, Prr[i + 1] - Prr[i]);
      ns += (Prr[i + 1] - Prr[i] > 1);
    }
    if (k > 1) {
      printf("   * %d strongly connected components, %d non-singleton, largest = %.1f %%\n", k, ns, 100.0 * largest / S->n);

      /* in place permute p to reflect block decomposition */
      for(i = 0; i < S->n; i++) {
	Sp[i] = p[ rr[1] + Sp[i] ];
      }
      for(i = 0; i < S->n; i++) {
	p[ rr[1] + i ] = Sp[i];
      }

      /* in place permute q to reflect block decomposition */
      for(j = 0; j < S->m; j++) {
	Sq[j] = q[ cc[2] + Sq[j] ];
      }
      for(j = 0; j < S->m; j++) {
	q[ cc[2] + j ] = Sq[j];
      }


    }

    spasm_csr_free(S);
    spasm_partition_free(P);
  }

  /* ------------------- V --------------------- */
  if (v) {
    V = spasm_submatrix(B, rr[2], rr[4], cc[3], cc[4], SPASM_IGNORE_VALUES);
    printf("*) V (%d x %d -- %d nnz) : \n",  V->n, V->m, spasm_nnz(V));

    Himatch = spasm_malloc(H->m * sizeof(int));
    // FIXME this is not finished
    for(j = cc[3]; j < cc[4]; j++) {
      Himatch[i] = qinv[ jmatch[ p[i] ] ];
    }
    P = spasm_connected_components(V, NULL, NULL, NULL);
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
      for(i = 0; i < V->n; i++) {
	Vp[i] = p[ rr[2] + Vp[i] ];
      }
      for(i = 0; i < V->n; i++) {
	p[ rr[2] + i] = Vp[i];
      }

      /* in place permute q to reflect block decomposition */
      for(j = 0; j < V->m; j++) {
	Vq[j] = q[ cc[3] + Vq[j] ];
      }
      for(j = 0; j < V->m; j++) {
	q[ cc[3] + j] = Vq[j];
      }

    }
    spasm_csr_free(V);
    spasm_partition_free(P);

    V = spasm_submatrix(B, rr[2], rr[3], cc[3], cc[4], SPASM_IGNORE_VALUES);
    P = spasm_strongly_connected_components(V);
    Prr = P->rr;
    k = P->nr;
    largest = -1;
    ns = 0;
    for(i = 0; i < k; i++) {
      largest = spasm_max(largest, Prr[i + 1] - Prr[i]);
      ns += (Prr[i + 1] - Prr[i] > 1);
    }
    if (k > 1) {
      printf("V' (%d x %d) : %d strongly connected components, %d non-singleton, largest = %.1f %%\n", V->n, V->m, k, ns, 100.0 * largest / V->n);
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
