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
 * type de données qui décrit un bloc rectangulaire.
 */
typedef struct {
  int i0; //  Le bloc est M[i0:i1, j0:j1], c'est-à-dire
  int j0; //     les lignes   de l'intervalle [i0; i1[
  int i1; //  et les colonnes de l'intervalle [j0; j1[
  int j1;
  int r;  // le rang du bloc
} block_t;


/*
 * type de données qui décrit un intervalle.
 */
typedef struct {
  int a; // décrit l'intervalle 
  int b; // [a; b[
} interval_t;


/*
 * type de données qui pointe vers un intervalle de ligne 
 * et un intervalle de colonne
 */
typedef struct {
  int i; // Pointe sur le i-ième intervalle de ligne
  int j; // et le j-ième intervalle de colonne.
} index_t;


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
 * Si blocks != NULL, stocke les blocs dans "blocks" à partir de l'indice *start (la valeur est modifiée)
 *
 * La position du coin inférieur droit du dernier bloc rencontrée doit être passée dans *last_i et *last_j. Ces valeurs sont modifiées.
 *
 */
void count_blocks(const spasm *M, spasm_cc *Y, block_t *blocks, int *start, int *last_i, int *last_j) {
  int i, j, a,b,c,d, r;

  for(i = 0; i < Y->CC->nr; i++) {
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	  a = Y->SCC[i]->rr[j];
	  b = Y->SCC[i]->cc[j];
	  c = Y->SCC[i]->rr[j + 1];
	  d = Y->SCC[i]->cc[j + 1];

	  if ((*last_i != a) || (*last_j != b)) {
	    // "hole" between last seen block and this one
	    if (blocks != NULL) {
	      blocks[*start].i0 = *last_i;
	      blocks[*start].j0 = *last_j;
	      blocks[*start].i1 = a;
	      blocks[*start].j1 = b;
	      blocks[*start].r  = 0;
	    }
	    (*start)++;
	  }

	  if (blocks != NULL) {
	    r = submatrix_rank(M, a, b, c, d);
	    blocks[*start].i0 = a;
	    blocks[*start].j0 = b;
	    blocks[*start].i1 = c;
	    blocks[*start].j1 = d;
	    blocks[*start].r  = r;
	  }
	  (*start)++;
	  (*last_i) = c;
	  (*last_j) = d;
      }
    }
  }
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
  int k, last_i, last_j;

  // étape 1 : détermine le nombre de blocs
  k = 0;
  last_i = 0;
  last_j = 0;
  if (DM->H != NULL) {
    count_blocks(M, DM->H, NULL, &k, &last_i, &last_j);
  }
  if (DM->S != NULL) {
    count_blocks(M, DM->S, NULL, &k, &last_i, &last_j);
  }
  if (DM->V != NULL) {
    count_blocks(M, DM->V, NULL, &k, &last_i, &last_j);
  }

  // étape 2 : allouer la liste des blocs
  *blocks = spasm_malloc(sizeof(block_t) * k);

  // étape 3 : remplir la liste des blocs
  k = 0;
  last_i = 0;
  last_j = 0;
  if (DM->H != NULL) {
    count_blocks(M, DM->H, *blocks, &k, &last_i, &last_j);
  }
  if (DM->S != NULL) {
    count_blocks(M, DM->S, *blocks, &k, &last_i, &last_j);
  }
  if (DM->V != NULL) {
    count_blocks(M, DM->V, *blocks, &k, &last_i, &last_j);
  }
  return k;
}


/*
 * Etant donnée le rang des blocks diagonaux, Bi, d'une matrice M, triangulaire par
 * blocs, détermine la liste des intervalles définis par [Bi.i0 + Bi.r ; Bi.i1[
 * d'une part (intervalles de lignes), et [Bi.j0 +Bi.r ; Bi.j1[ d'autre part
 * (intervalles de colonnes). 
 */
void intervals_fill(interval_t *R, interval_t *C, block_t *blocks, int n_blocks ) {
  int i;
  for (i=0; i<n_blocks; i++) {

    // Intervalles des lignes :
    R[i].a = blocks[i].i0 + blocks[i].r;
    R[i].b = blocks[i].i1;

    // Intervalles des colonnes :
    C[i].a = blocks[i].j0 +blocks[i].r;
    C[i].b = blocks[i].j1;
  }

}


void intervals_list(interval_t **R, interval_t **C, block_t *blocks, int n_blocks) {
  *R = spasm_malloc(sizeof(interval_t) * n_blocks);
  *C = spasm_malloc(sizeof(interval_t) * n_blocks);

  intervals_fill(*R, *C, blocks, n_blocks);
}


/*
 * Etant donnée les listes des intervalles "intéressants" de lignes et de colonnes, 
 * détermine le nombre de blocs à traîter sur la première diagonale supérieure de
 * M et remplit ces blocs.
 */

int count_other_blocks(const spasm *M, const interval_t *R, const interval_t *C, block_t *other_b, index_t *where, int n_blocks) {
  int k, i, j, nbl;
  nbl=0;
  for (k=n_blocks-1; k>0; k--) {
    if (C[k].a < C[k].b) {
      i = k-1;
      while ((R[i].a == R[i].b )&& (i>=0)) {
	i--;
      }
      if (i==-1) {
	return nbl;
      }
      j=i+1;
      while (C[j].a == C[j].b) {
	j++;
      }
      //remplissage des blocks.
      other_b[nbl].i0 = R[i].a;
      other_b[nbl].i1 = R[i].b;
      other_b[nbl].j0 = C[j].a;
      other_b[nbl].j1 = C[j].b;
      other_b[nbl].r = submatrix_rank(M, R[i].a, C[j].a, R[i].b, C[j].b);

      //Positions des blocks sur la diagonale.
      where[nbl].i = i;
      where[nbl].j = j;

      nbl++;
    }
  }
  return nbl;
}

int other_blocks_list(const spasm *M, const interval_t *R, const interval_t *C, block_t **other_b, index_t **where, int n_blocks) {
  int nbl, mem;
  mem = n_blocks-1;

  //allouer la mémoire des blocks grace au pointeurs temporaires;
  *other_b = spasm_calloc(mem, sizeof(block_t));
  *where = spasm_calloc(mem, sizeof(index_t));

  //remplir les blocks
  nbl = count_other_blocks(M, R, C, *other_b, *where, n_blocks);



  return nbl;
}


/*
 * Etant donnés les blocks d'une diagonale supérieure de M, détermine les intervalles intéressants
 * de la diagonale supérieure suivante.
 */
void change_intervals(interval_t *R, interval_t *C, const block_t *other_b, const index_t *where, int nbl) {
  int k;
  for (k=0; k<nbl; k++) {
    R[where[k].i].a = R[where[k].i].a + other_b[k].r;
    C[where[k].j].a = C[where[k].j].a + other_b[k].r;
  }
}


/*
 * Détermine la liste des blocks sur la k-ième diagonale supérieure de M
 */
int blocks_list_on_diag_k(const spasm *M, interval_t *R, interval_t *C, block_t *other_b, index_t *where, int n_blocks, int k) {
  int i, nbl;
  nbl = n_blocks;
  for (i=1; i<=k; i++) {
    nbl = other_blocks_list(M, R, C, &other_b, &where, nbl);
    change_intervals(R, C, other_b, where, nbl);
    free(other_b);
    free(where);
  }
  return nbl;
}


int main() {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n_blocks, i, *qinv, nbl, nbl1;
    block_t *blocks, *other_b;
    interval_t *R, *C;
    index_t *where;

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

    // Donne liste des intervalles de lignes et de colonnes.
    intervals_list(&R, &C, blocks, n_blocks);

    // Compte le nombre de blocs sur la première diagonale supérieure.
    nbl = other_blocks_list(B, R, C, &other_b, &where, n_blocks);

    // Compte le nombre de blocs sur la kième diagonale supérieure.
    nbl1 = blocks_list_on_diag_k(B, R, C, other_b, where, n_blocks, 1);

    // affichage
    for(i = 0; i < n_blocks; i++) {
      printf("%d ; %d ; % d; % d; %d \n", blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1, blocks[i].r);
      }
    printf("%d; %d; %d \n",n_blocks, nbl, nbl1);

    return 0;
}
