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
 * type de donnée qui décrit les blocs se présent sur une diagonale.
 */
typedef struct {
  int nbmax; // nombre de blocs sur la diagonale au maximum.
  int *nr; // numéro de l'intervalle ligne correspondant.
  int *nc; //numéro de l'intervalle colonne correspondant.
} diag_t;


/* alloue la mémoire d'un diag_t
 */
diag_t * diag_alloc(int nbmax) {
  diag_t *D;

  D = spasm_malloc(sizeof(diag_t));

  D->nbmax = nbmax;
  D->nr = spasm_calloc(nbmax, sizeof(int));
  D->nc = spasm_calloc(nbmax, sizeof(int));

  return D;
}

void diag_free(diag_t *D) {
  if(D == NULL) return;

  free(D->nr);
  free(D->nc);
  free(D);
}

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
  if (spasm_nnz(C) == 0) {
    spasm_csr_free(C);
    return 0;
  }

  // calcule la décomposition LU
  p = spasm_cheap_pivots(C);
  LU = spasm_LU(C, p, SPASM_DISCARD_L); // on se fiche de L
  free(p);

  // note le nombre de lignes non-nulles de U
  r = LU->U->n;

  // libère la sous-matrice et la mémoire dont on n'a plus besoin.
  spasm_free_LU(LU);
  spasm_csr_free(C);

  return r;
}


/*
 * Stocke le L de la décomposition LU dans L
 *  
 *
spasm * submatrix_L(const spasm *M, int a, int b, int c, int d) {
  spasm *C, *L;
  int *p;
  spasm_lu *LU;

  //extrait la sous_matrice
  C = spasm_submatrix(M, a, c, b, d, SPASM_WITH_NUMERICAL_VALUES);
  if (spasm_nnz(C)==0) {
    spasm_csr_free(C);
    L = NULL;
    return 0;
  }

  //calcule la décomposition LU
  p = spasm_cheap_pivots(C);
  LU = spasm_LU(C, p, SPASM_KEEP_L); // on garde L
  free(p);

  // récupère L
  L = LU->L;

  // libère la mémoire en trop.
  spasm_csr_free(C);
  spasm_csr_free(LU->U);
  free(LU->qinv);
  free(LU->p);

  return L;
  
}
*/


/*
 * Renvoie le rang de M[a:c, b:d].
 *
int submatrix_nnz(const spasm *M, int a, int b, int c, int d) {
  spasm *C;
  int r;

  // extrait la sous-matrice
  C = spasm_submatrix(M, a, c, b, d, SPASM_WITH_NUMERICAL_VALUES);
  r = spasm_nnz(C);
  spasm_csr_free(C);
  return r;
  }
*/


// incrémente empty de 1 si le bloc "block" est vide :
int is_block_empty(block_t block, int empty) {
  int k, r;
  r = block.r;
  k = ((r == 0) ? 1 : 0);
  empty = empty + k;
  return empty;
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
 *
 */


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

  // en principe on arrive forcément en bas de la matrice
  assert(last_i == M->n);

  // si on atteint pas le bord droit de la matrice, prévoir un bloc
  // vide pour aller jusqu'au bout
  if (last_j != M->m) {
    k++;
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

  // si on atteint pas le bord droit de la matrice, ajouter un bloc
  // vide pour aller jusqu'au bout
  if (last_j != M->m) {
    (*blocks)[k].i0 = last_i;
    (*blocks)[k].j0 = last_j;
    (*blocks)[k].i1 = M->n;
    (*blocks)[k].j1 = M->m;
    (*blocks)[k].r  = 0;
    k++;
  }

  return k;
}

/*
 * Etant donné un bloc diagonal, trouve le L de sa décomposition LU 
 * et calcule le produit du bas de l'inverse de L et du bloc "derrière"
 * le bloc diagonal  
 */
/*void find_salmon_block (const spasm *M, const block_t block, spasm *S) {
  spasm *L, *B;
  int i0, i1, j1, j2, from, to;

  i0 = block.i0;
  i1 = block.i1;
  j1 = block.j1;
  j2 = M->m;
  from = block.r;
  to = i1;

  // B sous-matrice définie associée au bloc [i0:i1, j1:j2]
  B = spasm_submatrix(M, i0, i1, j1, j2, SPASM_WITH_NUMERICAL_VALUES);
  L = submatrix_L(M, i0, j1, i1, j2);

  linvxm(L, B, from, to, S, SPASM_IDENTITY_PERMUTATION);

  }*/

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
 * Etant donnée une matrice diagonale par blocs, renvoie le
 * numéro du bloc diagonal correspondant à la colonne j pour tout j
 */
void column_diag_number(const spasm *M, const block_t *blocks, int *Q) {
  int j, m, k;

  m = M->m;
  k = 0;  

  for (j = 0; j < m; j++) {   
    while (blocks[k].j1 <= j) {
      k++;
    }
    Q[j] = k;
  }

}

/*
 * Etant donnés trois entiers l et c et une diagonale D, vérifie
 * si (l,c) est le dernier bloc ajouté dans D, si ce n'est pas le cas
 * ajouter (l,c) à D.
 * 
 * 
 */
int block_mark(diag_t *D, int l, int c, int last_b) {
  int *nr, *nc;

  nr = D->nr;
  nc = D->nc;

  if (last_b == -1) {
    last_b++;
    nr[last_b] = l;
    nc[last_b] = c;
  }
  else {
    if (nr[last_b] == l && nc[last_b] == c) return last_b; // <--- si bloc dernier marqué ne rien faire.
    last_b++;
    nr[last_b] = l;
    nc[last_b] = c;
  }
  return last_b;
}


/*
 * Pour chaque i dans l'intervalle de ligne l trouver j tel que
 * A[i,j] != 0, vérifier si le block d'intervalles (l, Q[j]) est
 * marqué comme appartenant à la diagonale Q[j] - l, et si ce n'est
 * pas le cas, le marquer.
 *
 * La valeur renvoyée correspond au nombre de blocs traités.
 */
int block_repartition(const spasm *M, const block_t *blocks, int n_block, const int *Q, diag_t **D, int *last_b) {
  int l, i, p, j, i0, i1, *Mp, *Mj, c, k, nbl;

  Mp = M->p;
  Mj = M->j;
  nbl = 0;

  for (l = 0; l < n_block; l++) {
    i0 = blocks[l].i0;
    i1 = blocks[l].i1;

    // trouver les emplacements des blocs non vides situés sur l'intervalle de lignes n°l.
    for (i = i0; i < i1; i++) {
      for (p = Mp[i]; p < Mp[i+1]; p++) {
	j = Mj[p];
	c = Q[j];
	k = c - l;
	last_b[k] = block_mark(D[k], l, c, last_b[k]);
      }
    }
  }

  // compter les blocs traités.
  for (k = 0; k < n_block; k++) {
    nbl = nbl + last_b[k];
  }
  return nbl;
}



/*
 * Etant donnée les listes des intervalles "intéressants" de lignes et de colonnes, 
 * détermine le nombre de blocs à traîter sur la première diagonale supérieure de
 * M et remplit ces blocs.
 */

int count_other_blocks(const spasm *M, const interval_t *R, const interval_t *C, block_t *other_b, const block_t *blocks, int n_blocks) {
  int k, i, j, nbl;
  nbl=0;
  for (k=0; k<n_blocks-1; k++) {
    if (R[k].a < R[k].b) {
     j = k+1;
      while ((C[j].a == C[j].b )&& (j<=n_blocks-1)) {
	j++;
      }
      if (j==n_blocks) {
	return nbl;
      }
      i=j-1;
      while (R[i].a == R[i].b) {
	i--;
      }
      //remplissage des blocks.
      other_b[nbl].i0 = blocks[i].i0;
      other_b[nbl].i1 = blocks[i].i1;
      other_b[nbl].j0 = blocks[j].j0;
      other_b[nbl].j1 = blocks[j].j1;
      other_b[nbl].r = submatrix_rank(M, blocks[i].i0, blocks[j].j0, blocks[i].i1, blocks[j].j1);

      nbl++;
    }
  }
  return nbl;
}

int other_blocks_list(const spasm *M, const interval_t *R, const interval_t *C, block_t **other_b, const block_t *blocks, int n_blocks) {
  int nbl, mem;
  mem = n_blocks-1;

  //allouer la mémoire des blocks grace au pointeurs temporaires;
  *other_b = spasm_calloc(mem, sizeof(block_t));

  //remplir les blocks
  nbl = count_other_blocks(M, R, C, *other_b, blocks, n_blocks);



  return nbl;
}

/*
 * Recopier une liste de blocks
 */
void blocks_dupli(block_t *blocks2, const block_t *blocks1, int nbl) {
  int i;
  for (i=0; i<nbl; i++) {
	blocks2[i].i0 = blocks1[i].i0;
	blocks2[i].i1 = blocks1[i].i1;
	blocks2[i].j0 = blocks1[i].j0;
	blocks2[i].j1 = blocks1[i].j1;
	blocks2[i].r = blocks1[i].r;
      }
}

void blocks_copy(block_t **blocks2, const block_t *blocks1, int nbl) {
  *blocks2 = spasm_malloc(nbl*sizeof(block_t));
  blocks_dupli(*blocks2, blocks1, nbl);
}


int main() {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n_blocks, i, j, *qinv, nbl, n_tot;
    block_t *blocks1, *blocks2;
    diag_t **D;
    //interval_t *R, *C;

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
    n_blocks = block_list(B, x, &blocks1);

    free(x);

    int k, nbmax, *last_b, *Q;
 
    // allocation mémoire de diag
    D = spasm_malloc(n_blocks * sizeof(diag_t));
  
    for (k = 0; k < n_blocks; k++) {
      nbmax = n_blocks - k;
      D[k] = diag_alloc(nbmax);
    }
   
    // allocation mémoire de last_b et Q
    last_b = malloc(n_blocks * sizeof(int));
    spasm_vector_set(last_b, 0, n_blocks, -1);
    
    Q = malloc(B->m * sizeof(int));

    // trouver le numéro de l'intervalle auquel appartient une colonne.
    column_diag_number(B, blocks1, Q);

   
    // trouver la répartition des blocs non vide.
    nbl = block_repartition(B, blocks1, n_blocks, Q, D, last_b); 

    //affichage
    n_tot = n_blocks + 1;
    n_tot = n_tot * n_blocks;
    n_tot = n_tot/2;

    printf("nombre total de blocs : %d\n", n_tot);
    printf("nombre de blocs non vides : %d\n", nbl);

    //libérer la mémoire
    free(blocks1);
    free(Q);
    free(last_b);
    for (k = 0; k < n_blocks; k++) {
      diag_free(D[k]);
    }
    free(D);
    spasm_csr_free(B);
    spasm_csr_free(A);


    //affichage
    //int nnz_diag = 0;
    //empty = 0;
    //for (i=0; i<n_blocks; i++) {
      //      printf("%d ; %d ; %d ; %d ; %d \n", blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1, blocks[i].r);
      //int dim_i =  blocks1[i].i1 - blocks1[i].i0;
      //int dim_j =  blocks1[i].j1 - blocks1[i].j0;
      //empty = is_block_empty(blocks1[i], empty);

	//printf("%d ; %d ; %d ; %d  ---> %d x %d    %d\n", blocks1[i].i0, blocks1[i].j0, blocks1[i].i1, blocks1[i].j1, dim_i, dim_j, r);
	//nnz_diag += r;
      
    //}

    // printf(" \n %d \n------------------------------\n", n_blocks);
    // printf("NNZ en tout : %d, NNZ sur la diagonale : %d\n", spasm_nnz(B), nnz_diag);
    // exit(0);

    // Compte le nombre de blocs sur les diagonales supérieures.
    //nbl = n_blocks;
    //n_tot = n_blocks;
    //for (j=0; j<10; j++) {

    // Détermine les intervalles de lignes et de colonnes.
    //intervals_list(&R, &C, blocks1, nbl);
    //blocks_copy(&blocks2, blocks1, nbl);
    //free(blocks1);

    // compte les blocks sur la diagonale.
    //nbl = other_blocks_list(B, R, C, &blocks1, blocks2, nbl);
    //n_tot = n_tot + nbl;

    // affichage
    /*
     *for (i=0; i<nbl; i++) {
	int dim_i =  blocks1[i].i1 - blocks1[i].i0;
	int dim_j =  blocks1[i].j1 - blocks1[i].j0;
	int r = blocks1[i].r;
	empty = is_block_empty(blocks1[i], empty);

	//printf("%d ; %d ; %d ; %d  ---> %d x %d    %d\n", blocks1[i].i0, blocks1[i].j0, blocks1[i].i1, blocks1[i].j1, dim_i, dim_j, r);
      }
    */
      //printf("%d \n------------------------------\n", empty);

    // libère la mémoire.
    //free (R);
    //free (C);
    //free (blocks2);


    //}
   
//printf("%d \n------------------------------\n", empty);

    // affichage
    return 0;
}
