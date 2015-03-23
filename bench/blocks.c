#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/*
 * Dans toute la suite, X[i:j] désigne les élements de X d'indice i <= ... < j,
 * C'est-à-dire de i inclus à j exclus. Il s'ensuit que i:j désigne l'intervalle [i; j[
 *
 * Cette notation bien pratique est héritée du langage Python.
 */


/*
 * type de données qui décrit un bloc carré.
 */
typedef struct {
  int i0; //  Le bloc est M[i0:i1, j0:j1], c'est-à-dire
  int j0; //     les lignes   de l'intervalle [i0; i1[
  int i1; //  et les colonnes de l'intervalle [j0; j1[
  int j1;
  int r;  // le rang du bloc
} block_t;


/*
 * Renvoie le rang de M[a:c, b:d].
 */
int submatrix_rank(const spasm *M, int a, int b, int c, int d) {
  spasm *C;
  int *p;
  spasm_lu *LU;
  int r;

  // extrait la sous-matrice
  C = spasm_submatrix(M, a, c, b, d, SPASM_WITH_NUMERICAL_VALUES);

  // calcule la décomposition LU
  p = spasm_cheap_pivots(C);
  LU = spasm_LU(C, p, SPASM_DISCARD_L); // on se fiche de L
  free(p);

  // note le nombre de lignes non-nulles de U
  r = LU->U->n;

  // libère la décomposition et la sous-matrice
  spasm_free_LU(LU);
  spasm_csr_free(C);

  return r;
}


/*
 * Renvoie "start" + (le nombre de blocs décrits par les composantes connexes dans "Y").
 *
 * Si blocks != NULL, stocke les blocs dans "blocks" à partir de l'indice start.
 *
 */
int count_blocks(const spasm *M, spasm_cc *Y, block_t *blocks, int start) {
  int i, j, a,b,c,d, r;

  for(i = 0; i < Y->CC->nr; i++) {
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	if (blocks != NULL) {
	  a = Y->SCC[i]->rr[j];
	  b = Y->SCC[i]->cc[j];
	  c = Y->SCC[i]->rr[j + 1];
	  d = Y->SCC[i]->cc[j + 1];
	  r = submatrix_rank(M, a, b, c, d);

	  blocks[start].i0 = a;
	  blocks[start].j0 = b;
	  blocks[start].i1 = c;
	  blocks[start].j1 = d;
	  blocks[start].r  = r;
	}
	start++;
      }
    }
  }
  return start;
}

/*
 * Etant donné une matrice M déjà permutée sous forme triangulaire par
 * blocs, et la decsription de la décomposition, détermine la liste
 * des blocs diagonaux.
 *
 * Ceci est rendu nécessaire par le fait que la structure de donnée
 * choisie pour stocker la décomposition n'est pas du tout pratique.
 *
 * le nombre de blocs est renvoyé. Il faut passer un pointeur vers une
 * liste de blocs, qui est modifiée.
 */
int block_list(const spasm *M, const spasm_dm *DM, block_t **blocks) {
  int k;

  // étape 1 : détermine le nombre de blocs
  k = 0;
  if (DM->H != NULL) {
    k = count_blocks(M, DM->H, NULL, k);
  }
  if (DM->S != NULL) {
    k = count_blocks(M, DM->S, NULL, k);
  }
  if (DM->V != NULL) {
    k = count_blocks(M, DM->V, NULL, k);
  }

  // étape 2 : allouer la liste des blocs
  *blocks = spasm_malloc(sizeof(block_t) * k);

  // étape 3 : remplir la liste des blocs
  k = 0;
  if (DM->H != NULL) {
    k = count_blocks(M, DM->H, *blocks, k);
  }
  if (DM->S != NULL) {
    k = count_blocks(M, DM->S, *blocks, k);
  }
  if (DM->V != NULL) {
    k = count_blocks(M, DM->V, *blocks, k);
  }
  return k;
}


int main() {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n_blocks, i, *qinv;
    block_t *blocks;

    // charge la matrice depuis l'entrée standard
    T = spasm_load_sms(stdin, 42013);
    A = spasm_compress(T);
    spasm_triplet_free(T);

    // calcule la décomposition de A
    x = spasm_dulmage_mendelsohn(A);

    // B = A permutée sous forme triangulaire par blocs
    qinv = spasm_pinv(x->DM->q, A->m);
    B = spasm_permute(A, x->DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
    free(qinv);

    // calcule la liste des blocs
    n_blocks = block_list(B, x, &blocks);

    // affichage
    for(i = 0; i < n_blocks; i++) {
      printf("%d ; %d ; % d; % d; %d \n", blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1, blocks[i].r);
    }
    return 0;
}
