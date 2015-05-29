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
 * Type de données qui décrit un intervalle.
 */
typedef struct {
  int a; // décrit l'intervalle 
  int b; // [a; b[
} interval_t;


/*
 * Type de donnée qui donne la position d'un block.
 */
typedef struct {
  int c; // numéro de la colonne auquel le bloc appartient.
  int r; // numéro de l'interval de lignes auquel le bloc appartient.
} blk_t;


/*
 * Type de donnée qui décrit une matrice triangulaire supérieure,
 * dont on se moque de la valeur numérique des entrées.
 */
typedef struct {
  int nzmax; // nombre maximum d'entrée non nulle.
  int n; // nombre de ligne = nombre de diagonale.
  int *d; // pointeurs sur les diagonales
  int *i; // positions sur les lignes.
  int *j; // positions sur les colonnes.
} uptri_t;


/*
 * Type de donnée qui décrit les blocs se présent sur une diagonale.
 */
typedef struct {
  int ndiag; // numero de la diagonale.
  int nbmax; // nombre de blocs sur la diagonale au maximum.
  int *nr; // numéro de l'intervalle ligne correspondant.
  int *nc; //numéro de l'intervalle colonne correspondant.
} diag_t;


/*
 * Type de donnée qui définie les actions prévues en une ligne/colonne.
 * Utilise les listes simplement chaînées
 */
typedef struct action 
{
  int act; // élément de la liste chaînée
  struct action *next; // pointeur sur l'élément suivant
} action_t; 


/*
 * Type de donnée qui définie les actions prévues en unz ligne colonne.
 * Utilise les arbres binaires de recherche.
 */
typedef struct actree {
  int act; // valeur du noeud.
  struct actree *right; // pointeur sur la valeur de droite (supérieure)
  struct actree *left; // pointeur sur la valeur de gauche (inférieure)
} actree_t;

/* alloue la mémoire d'un diag_t
 */
diag_t * diag_alloc(int nbmax) {
  diag_t *D;

  D = spasm_malloc(sizeof(diag_t));

  D->ndiag = 0;
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
 * aloue la mémoire d'un uptri_t;
 */
uptri_t * uptri_alloc(int nzmax, int n, int look_row, int look_col) {
  uptri_t *T;

  T = spasm_malloc(sizeof(uptri_t));

  T->nzmax = nzmax;
  T->n = n;
  T->d = spasm_malloc( (n+1) * sizeof(int));
  T->i = (look_row ? spasm_malloc(nzmax * sizeof(int)) : NULL);
  T->j = (look_col ? spasm_malloc(nzmax * sizeof(int)) : NULL);


  return T;
}

/*
 * libère la mémoire d'un uptri_t;
 */
void uptri_free(uptri_t *T) {
  if (T == NULL) return;

  free(T->d);
  free(T->i);
  free(T->j);
  free(T);
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


// incrémente fill de 1 si le bloc "block" n'est pas vide :
int is_block_empty(block_t block, int fill) {
  int k, r;
  r = block.r;
  k = ((r == 0) ? 0 : 1);
  fill = fill + k;
  return fill;
}



/*
 * Renvoie "start" + (le nombre de blocs décrits par les composantes connexes dans "Y").
 *
 * Si blocks != NULL, stocke les blocs dans "blocks" à partir de l'indice *start (la valeur est modifiée)
 *
 * La position du coin inférieur droit du dernier bloc rencontrée doit être passée dans *last_i et *last_j. Ces valeurs sont modifiées.
 *
 */
void count_blocks(spasm_cc *Y, block_t *blocks, int *start) {
  int i, j, a,b,c,d;

  for(i = 0; i < Y->CC->nr; i++) {
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	  a = Y->SCC[i]->rr[j];
	  b = Y->SCC[i]->cc[j];
	  c = Y->SCC[i]->rr[j + 1];
	  d = Y->SCC[i]->cc[j + 1];

	  if (blocks != NULL) {
	    blocks[*start].i0 = a;
	    blocks[*start].j0 = b;
	    blocks[*start].i1 = c;
	    blocks[*start].j1 = d;
	  }
	  (*start)++;
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
int block_list(const spasm *M, const spasm_dm *DM, block_t **blocks_ptr) {
  int i, k;
  block_t *blocks;
  
  // étape 1 : détermine le nombre de blocs
  k = 0;
  if (DM->H != NULL) {
    count_blocks(DM->H, NULL, &k);
  }
  if (DM->S != NULL) {
    count_blocks(DM->S, NULL, &k);
  }
  if (DM->V != NULL) {
    count_blocks(DM->V, NULL, &k);
  }

  // étape 2 : allouer la liste des blocs
  blocks = spasm_malloc(sizeof(block_t) * k);
  *blocks_ptr = blocks;
  
  // étape 3 : remplir la liste des blocs
  k = 0;
  if (DM->H != NULL) {
    count_blocks(DM->H, blocks, &k);
  }
  if (DM->S != NULL) {
    count_blocks(DM->S, blocks, &k);
  }
  if (DM->V != NULL) {
    count_blocks(DM->V, blocks, &k);
  }

  // étape 4 : "malaxer" la liste des blocs (il faut que leurs coins se touchent).
  blocks[0].i0 = 0;
  blocks[0].j0 = 0;
  for (i = 1; i < k; i++) {
    blocks[i - 1].i1 = blocks[i].i0;
    blocks[i - 1].j1 = blocks[i].j0;
  }
  blocks[k-1].i1 = M->n;
  blocks[k-1].j1 = M->m;

  // étape 5 : calculer les rangs
  for (i = 0; i < k; i++) {
    blocks[i].r = submatrix_rank(M, blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1);
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
 * et celui correspondant à la ligne i pour tout i;
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
 * vérifie si un entier k appartient à la liste D et si ce n'est pas le cas l'ajouter
 *
 * la valeur renvoyée est le nombre d'éléments dans la liste.
 *
int diag_list(int *D, int k, int size) {
  int i, j, start, end, found;

  // initialisation de la liste.
  if (size == 0) {
    D[0] = k;
    return 1;
  }
  
  // vérifier si k appartient à la liste.
  
  //par dichotomie :

  start = 0;
  end = size;

  while( (end - start) > 1) {
    i = (start + end)/2;
    if (D[i] == k) return size; // si k est dans la liste ne rien faire.

    if (D[i] > k) end = i;
    else start = i;
  }

  if (D[start] == k) return size; // si k dans la liste ne rien faire.

  // Si k n'est pas dans la liste l'insérer.
  j = size;
  while (j > 0 && k < D[j - 1]) {
    D[j] = D[j-1];
    j--;
  }
  D[j] = k;
  size++;

  return size;
} 

 *
 * donne le nombre de diagonale qui contiennent au moins un block non vide.
 *
int diag_count(const spasm *M, const block_t *blocks, const int *Q, int *D) {
  int l, i, p, j, *Mp, *Mj, n, c, k, n_diag;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  l = 0;
  n_diag = 0;

    for (i = 0; i < n; i++) {
      while (blocks[l].i1 <= i) l++;

      for (p = Mp[i]; p < Mp[i+1]; p++) {
	j = Mj[p];
	c = Q[j];
	k = c - l;
	n_diag = diag_list(D, k, n_diag);
	printf("---%d---\n", n_diag);
      }
    }

  return n_diag;
}
*/

/*
 * Détermine le nombre de blocs non vide et les diagonales auquels ils appartiennent.
 *
 * Ici seul le dernier bloc observé à la diagonale k est stocké à D[k]. On n'a pas
 * d'information sur les autres blocs.
 */
int count_filled_blocks(const spasm *M, const block_t *blocks, int n_blocks, const int *Q) {
 int i, l, p, j, *Mp, *Mj, n, c, k, fill, *D;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  l = 0; // <--- numéro du bloc de ligne qu'on regarde.
  fill = 0;

  D = spasm_malloc(n_blocks * sizeof(int));
  for (i = 0; i < n_blocks; i++) {
    D[i] = -1;
  }

  for (i = 0; i < n; i++) {
    while (blocks[l].i1 <= i) l++;

    for (p = Mp[i]; p < Mp[i+1]; p++) {
      j = Mj[p];
      c = Q[j]; 
      k = c - l; // <--- numéro de la diagonale à laquelle appartient l'entrée.
      if (D[k] != l) {
	fill++;
	D[k] = l; // <--- l remplace la valeur précédente de D[k].
      }
    }
  }
  free(D);
  return fill;
}

/*
 * Détermine la liste des emplacements des blocs non vides si ils sont "intéressants"
 * la valeur renoyée est le nombre de blocs "intéressants" à regarder
 */
int filled_blocks_list(const spasm *M, const block_t *blocks, int n_blocks, const int *Q, blk_t *where) {
  int i, l, p, j, *Mp, *Mj, n, c, k, fill, *D, count;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  l = 0; // <--- numéro du bloc de ligne qu'on regarde.
  fill = 0;
  count = 0;

  D = spasm_malloc(n_blocks * sizeof(int));
  for (i = 0; i < n_blocks; i++) {
    D[i] = -1;
  }

  for (i = 0; i < n; i++) {
    while (blocks[l].i1 <= i) l++;

    for (p = Mp[i]; p < Mp[i+1]; p++) {
      j = Mj[p];
      c = Q[j]; 
      k = c - l; // <--- numéro de la diagonale à laquelle appartient l'entrée.
      if (D[k] != l) {
	fill++;
	D[k] = l; // <--- l remplace la valeur précédente de D[k].
	if ((k == 0) || (blocks[c].r < blocks[c].j1 - blocks[c].j0)) {
	  where[count].c = c;
	  where[count].r = l;
	  count++;
	 }
      }
    }
  }
  free(D);
  return count;
}


/*
 * Donne la matrice des positions des blocs sous forme uptri_t en regardant les lignes
 * à partir de la liste "where", du nombre de blocs "count" et
 * du nombre de blocs diagonaux "n_blocks".
 */
uptri_t * position_matrix_row_view(const blk_t *w, int n_blocks, int count) {
  int k, *wd, *Bd, *Bi, *tmp, sum, p;
  uptri_t *B;

  /* Pour chaque bloc w[k] trouver sa diagonale. */

  wd = spasm_malloc(count * sizeof(int)); // Allocation mémoire diagonale

  for (k = 0; k < count; k++) {
    wd[k] = w[k].c - w[k].r; // <--- diagonale à laquelle appartient le bloc w[k]
  }

  /* allocation du résultat */
  B = uptri_alloc(count, n_blocks, 1, 0); // <--- On regarde les lignes

  tmp = spasm_calloc(n_blocks, sizeof(int));

  Bd = B->d;
  Bi = B->i;

  /* compte le nombre d'entrée sur une diagonale */
  for (k = 0; k < count; k++) {
    tmp[ wd[k] ]++;
  }

  /* implémente les pointeurs de diagonales */
  sum = 0;
  for (k = 0; k < n_blocks; k++) {
    Bd[k] = sum;
    sum += tmp[k];
    tmp[k] = Bd[k];
  }
  Bd[n_blocks] = sum;

  /* trouve l'emplacement des entrées sur chaque diagonale */
  for (k = 0; k < count; k++) {
    p = tmp[ wd[k] ]++; // <-- p-ième entrée de B, sur la diagonale wd[k].
    Bi[p] = w[k].r;  // On note le numéro de ligne correspondant à l'entrée.
 }

  /* libération de mémoire, retourner le resultat. */
  free(tmp);
  free(wd);
  return B;
}

/*
 * Donne la matrice des positions des blocs sous forme uptri_t en regardant les colonnes
 * à partir de la liste "w", du nombre de blocs "count" et du nombre de blocs diagonaux
 * "n_blocks".
 *
 * Les colonnes sont triées par ordres croissant.
 */
uptri_t * position_matrix_col_view(const blk_t *w, int n_blocks, int count) {
  int k, *wd, *Bd, *Bj, *tmp, sum, p, tmpj, l;
  uptri_t *B;

  /* Pour chaque bloc w[k] trouver sa diagonale. */
  wd = spasm_malloc(count * sizeof(int));

  for(k = 0; k < count; k++) {
    wd[k] = w[k].c - w[k].r; 
  }

  /* allocation du résultat */
  B = uptri_alloc(count, n_blocks, 0, 1); // <--- On regarde les colonnes.

  tmp = spasm_calloc(n_blocks, sizeof(int));

  Bd = B->d;
  Bj = B->j;

  /* compte le nombre d'entrée sur une diagonale */
  for (k = 0; k < count; k++) {
    tmp[ wd[k] ]++;
  }

  /* implémente pointeurs de diagonales. */
  sum = 0;
  for(k = 0; k < n_blocks; k++) {
    Bd[k] = sum;
    sum += tmp[k];
    tmp[k] = Bd[k]; 
  }
  Bd[n_blocks] = sum;

  /* trouve les emplacements des entrées sur chaque diagonales (numéro de la colonne) */
  for (k = 0; k < count; k++) {
    p = tmp[ wd[k] ]++;
    Bj[p] = w[k].c;
  }

  /* trie les emplacements par ordre croissant (tri par insertion) */
  for (k = 0 ; k < n_blocks; k++) {
    for (p = Bd[k] + 1; p < Bd[k+1]; p++) {
      tmpj = Bj[p];
      l = p-1;
      while(l >= 0 && Bj[l] > tmpj) {
	Bj[l+1] = Bj[l];
	l--;
      }
      Bj[l+1] = tmpj;
    }
  }

  /* libération de mémoire auxiliaire, renvoi du résultat. */
  free(tmp);
  free(wd);

  return B;
}


/*
 * Pour un bloc donné par son entrée k et sa diagonale dans la matrice, 
 * si le bloc est susceptible d'être éliminé, alors, 
 * Er[r] = k (on regarde les lignes) et Ec[c] = k (on regarde les colonnes).
 * renvoie 1 si le bloc est succeptible d'être éliminé, 0 sinon.
 */
int may_be_eliminated(int k, int diag, const uptri_t *B, const block_t *blocks, int *Er, int *Ec) {
  int r, c, i0, i1, *Bi;

  Bi = B->i;

  r = Bi[k];
  c = Bi[k] + diag;

  i0 = blocks[r].i0;
  i1 = blocks[r].i1;

  if (i0 < i1) {
    Er[r] = k;
    Ec[c] = k;
    return 1;
  }
  return 0;
}

/*
 * Pour un bloc donné par son entrée k dans la matrice, 
 * trouve le dernier bloc éliminé à gauche.
 */
int left_elimination(int k, const int *Er, const uptri_t *B) {
  int *Bi;

  Bi = B->i;
 

  return Er[Bi[k]];

}

/*
 * Pour un bloc donné par son entrée k et sa diagonale diag, 
 * trouve le dernier bloc éliminé en dessous
 */
int under_elimination(int k, int diag, const int *Ec, const uptri_t *B) {
  int  *Bi, c;

  Bi = B->i;
  c = Bi[k] + diag;

  return Ec[c];
}


/*
 * Pour un indice k correspondant à une entrée sur la matrice dans une diagonale inférieure à diag
 * renvoie le numéro de la diagonale sur laquelle se trouve le bloc d'indice k.
 */
int diag_index(const uptri_t *B, int k, int diag) {
  int d, *Bd;

  Bd = B->d;

  assert(Bd[diag] > k);

  for (d = diag; d >=0 && Bd[d] > k; d--); // <--- numéro de la diagonale sur laquelle est le bloc d'indice k.

  return d ;
}


/*
 * Etant donnés r, et c, regarde si le bloc (r,c) correspond à une entrée de la matrice des positions.
 * renvoie si 0 si c'est le cas et 1 sinon.
 */
int is_entry_zero(const uptri_t *B, int r, int c) {
  int d, *Bi, *Bd, k;

  Bi = B->i;
  Bd = B->d;

  d = c - r;

  for (k = Bd[d]; k < Bd[d+1] && Bi[k] != r; k++);

  if (Bi[k] != r) return 1;

  return 0;
}


/*
 * Retourne 0 si k appartient à la liste chaînée list et 1 sinon.
 */
int is_action_new(action_t *list, int k) {
  action_t *tmp;

  tmp = list;

  while(tmp != NULL) {
    if(tmp->act == k) return 0;
    tmp = tmp->next;
  }

  return 1;
}


/*
 * Ajoute une action en tête de liste. 
 */
action_t * new_action(action_t *list, int k) {
  action_t *new;

  new = spasm_malloc(sizeof(action_t)); //<--- création d'un nouvel élément de la chaine
  new->act = k; // <--- la valeur de l'action est k.
  new->next = list; // <--- l'élément suivant est le premier élément de la liste list.

  return new;
}

/*
 * Etant donné c, parcourt les actions r d'une liste chaînée list, 
 * et regarde pour tout r, si (r,c) est dans la matrice B. 
 * Renvoie le total des actions déclanchées.
 */
int row_action_set_off(const uptri_t *B, action_t *list, int c) {
  int r, tot;
  action_t *tmp;

  tmp = list;
  tot = 0;

  while(tmp != NULL) {
    r = tmp->act;
    tot += is_entry_zero(B, r, c);
    tmp = tmp->next;
  }

  return tot;
}

/*
 * Etant donné r, parcourt les actions c d'une liste chaînée et regarde pour tout c
 * si (r, c) est dans la matrice B. Renvoie le total des actions déclanchées.
 */
int col_action_set_off(const uptri_t *B, action_t *list, int r) {
  int c, tot;
  action_t *tmp;

  tmp = list;
  tot = 0;

  while(tmp != NULL) {
    c = tmp->act;
    tot += is_entry_zero(B, r, c);
    tmp = tmp->next;
  }

  return tot;
}

/*
 * vide une liste d'action et libère la mémoire
 */
action_t * clear_action(action_t *list) {
 
  if (!list) return NULL;

  action_t *tmp;

  tmp = list->next;
  free(list);
  return clear_action(tmp);

}

/*
 * Ajoute une action dans un arbre binaire tree si celle-ci n'existe pas déjà.
 * renvoie 1 si l'action a été ajouté et 0 sinon.
 */
int new_act(actree_t **tree, int k) {
  actree_t *new, *tmp1, *tmp2;

  //Allocation et initialisation de l'arbre new.
  new = malloc(sizeof(actree_t));

  new->act = k;
  new->left = NULL;
  new->right = NULL;

  tmp1 = *tree;

  if (tmp1 == NULL) {
    *tree = new;
    return 1; // ajout du premier élément de la liste.
  }

  while(tmp1 != NULL) {

    if (k == tmp1->act) { 
      free(new);      
      return 0;
    }

    tmp2 = tmp1;
    if (k > tmp1->act) {
      tmp1 = tmp1->right;
      if(tmp1 == NULL) {
	tmp2->right = new;
	return 1;
      }
    }
    else {
      tmp1 = tmp1->left;
      if(tmp1 == NULL) {
	tmp2->left = new;
	return 1;
      }
    }
  }

  printf("Error : tree search fail for k = %d \n", k);
  return -1;
}

/*
 * Etant donné c, parcourt les actions r d'un arbre binaire tree, et regarde pour tout r, 
 * si (r,c) est dans la matrice B. Renvoie le total des actions déclanchées.
 */
void row_act_set_off(const uptri_t *B, actree_t *tree, int c, int *tot) {
  int r;
  actree_t *tmp;

  tmp = tree;

  if(!tmp) return;
  if(tmp->left) row_act_set_off(B, tmp->left, c, tot);
  if(tmp->right) row_act_set_off(B, tmp->right, c, tot);

  r = tmp->act;
  (*tot) += is_entry_zero(B, r, c);

}

/*
 * Etant donné r, parcourt les actions c d'un arbre binaire tree, et regarde pour tout c, 
 * si (r,c) est dans la matrice B. Renvoie le total des actions déclanchées.
 */
void col_act_set_off(const uptri_t *B, actree_t *tree, int r, int *tot) {
  int c;
  actree_t *tmp;
  
  tmp = tree;

  if (!tmp) return;
  if (tmp->left) col_act_set_off(B, tmp->left, r, tot);
  if (tmp->right) col_act_set_off(B, tmp->right, r, tot);

  c = tmp->act;
  (*tot) += is_entry_zero(B, r, c);
}

/*
 * Vide un arbre d'actions, libère la mémoire.
 */
void clear_act(actree_t **tree) {
  actree_t *tmp;

  if(*tree == NULL) return;

  tmp = *tree;
  if(tmp->left != NULL) clear_act(&tmp->left);
  if(tmp->right != NULL) clear_act(&tmp->right);

  free(tmp);
  *tree = NULL;

}


/*
 * Parcourt la matrice diagonale par diagonale, pour chaque entrée, stocke le bloc éliminé 
 * immediatement à gauche et celui immédiatement à droite. (On regarde les entrée à partir
 * de la première diagonale supérieure)
 *
 * Détermine si le bloc est susceptible d'être éliminé ou non
 * 
 * Ajoute les actions aux bonnes colonnes et aux bonnes lignes
 *
 * Effectue les actions si nécessaire.
 *
 * renvoie la nouvelle valeur du nombre de bloc
 */
int emergence_simulation(uptri_t *B, const block_t *blocks, int n_blocks, int *c_act, int *r_act) {
  int k, i, j, d, diag, count, *Bd, *Bi, *Ec, *Er,  *left, *under, l, u, elim, add, start; 
    //*c_act, *r_act;
  blk_t blk;
  // action_t **Row, **Col;
  actree_t **Rtree, **Ctree;

  Bd = B->d;
  Bi = B->i;
  count = B->nzmax;
  start = 0; // nombre d'action déclenchées, initialisée à 0 au début du programme.
 
  /* Allocation des espaces de mémoire.
   */

  // Allocation mémoire des tableaux de listes chaînées.
  //Row = spasm_malloc(n_blocks * sizeof(action_t*));
  //Col = spasm_malloc(n_blocks * sizeof(action_t*));

  //Allocation mémoire des tableaux d'arbres binaires.
  Rtree = spasm_malloc(n_blocks * sizeof(actree_t*));
  Ctree = spasm_malloc(n_blocks * sizeof(actree_t*));

  // Initialisation des listes chaînées :
  for (k = 0; k <n_blocks; k++) {
    //Row[k] = NULL; // Au départ une liste chaînées pointe sur NULL
    //Col[k] = NULL;

    Rtree[k] = NULL; // Au départ, les arbres pointent sur NULL
    Ctree[k] = NULL; 
  }


  // Allocation mémoire pour listes des blocs ayant subit une élimination.
  Ec = spasm_malloc(n_blocks * sizeof(int));
  Er = spasm_malloc(n_blocks * sizeof(int));

  // Allocation mémoire de under et left.
  left = spasm_malloc(count * sizeof(blk_t)); // pour tout k, left[k] désigne l'indice du dernier bloc éliminé à gauche du bloc d'indice k.
  under = spasm_malloc(count *sizeof(blk_t)); // under[k] désigne l'indice du dernier bloc éliminé sous le bloc d'indice k

  // Allocation mémoire de c_act et r_act, comptent les actions sur les colonnes et les lignes.
  // r_act = spasm_malloc(n_blocks * sizeof(blk_t)); // pour tout i, r_act[i] désigne le nombre d'action à effectuer à la ligne i.
  // c_act = spasm_malloc(n_blocks * sizeof(blk_t)); // pour tout j, c_act[j] désigne le nombre d'action à effectuer à la colonne j.


  /* Initialisation des données
   */
  
  for(k = 0; k < n_blocks; k++) {
    left[k] = -1; // Si k désigne un bloc sur la diagonale principale, il n'y a pas de bloc éliminé avant
    under[k] = -1; // on initialise à -1.
    
    Er[Bi[k]] = k; // Pour chaque ligne, on élimine l'entrée correspondante
    Ec[Bi[k]] = k; // Pour chaque colonne, on élimine l'entrée correspondante

    r_act[k] = 0; // Au début du programme, aucune action n'est prévue.
    c_act[k] = 0;

  }


  /* Parcours de la matrice pour les diagonales supérieures.
   */

  for (diag = 1; diag < n_blocks; diag++) {
   
    for (k = Bd[diag]; k < Bd[diag + 1]; k++) {
      blk.r = Bi[k];
      blk.c = Bi[k] + diag;

      /* Donne l'entrée correspondant au dernier bloc éliminé à gauche
       * Et en dessous du bloc (i, j)
       */
      
      left[k] = left_elimination(k, Er, B);
      under[k] = under_elimination(k, diag, Ec, B);

      // if (left[k] == -1 || under[k] == -1) {
      //printf("Error entry %d of position matrix \n", k);
      //return 0;
      //}
      
      /* Regarde les actions à ajouter 
       * sur les colonnes à gauche.
       */
      l = left[k];
      d = diag;
      while (l != -1) {

	d = diag_index(B, l, d); // détermine les diagonales des bloc éliminés à gauche du bloc d'indice k.
	j = d + Bi[l]; // colonne du bloc d'indice l

	// teste si blk.c appartient à la liste des actions prévues en j
	//add = is_action_new(Col[j], blk.c);

	// ajoute blk.c dans la liste des actions prévues en j, si ce n'est pas le cas.
	//if(add) {    
	  //Col[j] = new_action(Col[j], blk.c);
	//	}

	// incrémente le compteur si nécessaire.
	//	c_act[j] += add;

	c_act[j] += new_act(&Ctree[j], blk.c);

	l = left[l];
      } 

      //On compte les actions que le bloc déclenche sur la ligne blk.r
      //start += row_action_set_off(B, Row[blk.r], blk.c);
       row_act_set_off(B, Rtree[blk.r], blk.c, &start);

      /* Regarde si le bloc doit être éliminé ou non.
       * met à jour les listes Er et Ec
       */
      elim = may_be_eliminated(k, diag, B, blocks, Er, Ec);

      /* Si le bloc doit être éliminé on regarde les action à 
       * ajouter sur les lignes en dessous.
       */
      if(elim == 1) {
	u = under[k];
	while(u != -1) {
	  i = Bi[u]; // ligne correspondant au bloc d'indice u

	  // teste si blk.r appartient à la listes des actions prévues en i.
	  //add = is_action_new(Row[i], blk.r);

	  // ajoute blk.r dans la liste d'action prévues en i, si ce n'est pas le cas.
	  //if(add) {
	     //Row[i] = new_action(Row[i], blk.r);
	  //	  }
	  
	  // incrémente le compteur d'action de la ligne i.
	  //r_act[i] += add;

	  r_act[i] += new_act(&Rtree[i], blk.r);
	  
	  u = under[u];
	}

	//On compte les actions que le bloc déclenche sur la colonne blk.c
	//start += col_action_set_off(B, Col[blk.c], blk.r);
	col_act_set_off(B, Ctree[blk.c], blk.r, &start);
      }
      
      /* Regarde si le bloc blk ne déclenche pas lui-même d'action.
       * Dans tous les cas, on déclenche les actions prévues sur la ligne blk.r
       * Si on élimine le bloc, on déclenche les actions prévues sur la colonne blk.c
       */

    }
    
  }

  /* Libération de la mémoire auxiliaire.
   */
   
  for(k = 0; k < n_blocks; k++) {
    //Row[k] = clear_action(Row[k]);
     //Col[k] = clear_action(Col[k]);
    
    clear_act(&Rtree[k]);
    clear_act(&Ctree[k]);
  }
  
  
  free(Ec);
  free(Er);
  free(left);
  free(under);
  //free(Row);
  //free(Col);
  free(Rtree);
  free(Ctree);
  //free(r_act);
  //free(c_act);

  return start;

}


/*
 * Etant donnés deux entiers l et c et une diagonale D, vérifie
 * si (l,c) est le dernier bloc ajouté dans D, si ce n'est pas le cas
 * ajouter (l,c) à D.
 * 
 * 
 *
int block_mark(diag_t *diags, int l, int c, int k, int next_b) {
  int *nr, *nc;

  nr = diags->nr;
  nc = diags->nc;

  diags->ndiag = k;

  if (next_b == 0) {
    nr[next_b] = l;
    nc[next_b] = c;
    next_b++;
   return  next_b;
  }
   
  if (nr[next_b - 1] == l && nc[next_b - 1] == c) return next_b; // <--- si bloc dernier marqué ne rien faire.

  nr[next_b] = l;
  nc[next_b] = c;
  next_b++;
  return next_b;
}

 *
 * Etant donnée une liste L de taille t, d'entiers disincts et 
 * ordonnés, et un élément k de L, trouve l'unique i vérifiant
 * L[i] = k.
 *
int dicho(const int *L, int k, int t) {
  int i, start, end, found;

  found = 0; // <--- On n'a pas encore trouvé i. 
  start = 0;
  end = t;
    
  while (found == 0 && (end - start) > 1) {
    i = (start + end)/2 ; 
    found = ((L[i] == k) ? 1 : 0); // <--- test si i est l'indice qu'on recherche.
    if (L[i] > k) end = i;
    else start = i;
  }

  if (L[start] == k) return start;
  else return -1; // valeur renvoyée par défaut si k n'appartient pas à la liste L.

}

 *
 * Pour chaque i dans l'intervalle de ligne l trouver j tel que
 * A[i,j] != 0, vérifier si le block d'intervalles (l, Q[j]) est
 * marqué comme appartenant à la diagonale Q[j] - l, et si ce n'est
 * pas le cas, le marquer.
 *
 * La valeur renvoyée correspond au nombre de blocs traités.
 *
int blocks_repartition(const spasm *M, const block_t *blocks, const int *Q, diag_t **diags, int n_diag, int *next_b, const int *D) {
  int l, i, p, j, n, *Mp, *Mj, c, k, nbl, d;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  nbl = 0;
  l = 0;
 

    // trouver les emplacements des blocs non vides situés sur l'intervalle de lignes n°l.
    for (i = 0; i < n; i++) {
      while (blocks[l].i1 <= i) l++;
      for (p = Mp[i]; p < Mp[i+1]; p++) {
	j = Mj[p];
	c = Q[j];
	k = c - l;

	// trouver l'unique d tel que D[d] = c - l.
	d = dicho(D, k, n_diag);
	if (d == -1) {
	  printf ("Error : digonal %d is supposed to be empty\n", k);
	  return -1;
	}
	next_b[d] = block_mark(diags[d], l, c, k, next_b[d]);
       
      }
    }


  // compter les blocs traités.
  for (k = 0; k < n_diag; k++) {
    nbl = nbl + next_b[k];
  }
  return nbl;
}

*/

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
    int n_blocks, i, *qinv;
    block_t *blocks1;
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
    for(i = 0; i < n_blocks; i++) {
      printf("%d : (%d, %d) -- (%d, %d)\n", i, blocks1[i].i0, blocks1[i].j0, blocks1[i].i1, blocks1[i].j1);
    }
    printf("blocs diagonaux : %d\n", n_blocks);

    free(x);

    int *Q, fill, count, *c_act, *r_act, start;
    blk_t *where;
    uptri_t *P, *S;
 

    // allocation mémoire de not_empty et Q
    //not_empty = malloc(n_blocks * sizeof(int));
    Q = malloc(B->m * sizeof(int));

    // trouver le numéro de l'intervalle auquel appartient une colonne.
    column_diag_number(B, blocks1, Q);
    printf("-------------------------\n");

    // trouver le nombre de diagonales ayant au moins un bloc non vide.
    // n_diags = diag_count(B, blocks1, Q, not_empty);

    //printf("%d\n", n_diags);


    fill = count_filled_blocks(B, blocks1, n_blocks, Q); // <--- nombre total de blocs non vide.

    printf("%d\n", fill);

    //allocation de mémoire de where.
    where = malloc(fill * sizeof(blk_t));
    for (i = 0; i < fill; i++) {
      where[i].c = -1;
      where[i].r = -1;
    }

    count = filled_blocks_list(B, blocks1, n_blocks, Q, where);
    P = position_matrix_row_view(where, n_blocks, count);
    S = position_matrix_col_view(where, n_blocks, count);

    printf("%d\n", count);


   
    r_act = spasm_malloc(n_blocks * sizeof(blk_t)); // pour tout i, r_act[i] désigne le nombre d'action à effectuer à la ligne i.
    c_act = spasm_malloc(n_blocks * sizeof(blk_t)); // pour tout j, c_act[j] désigne le nombre d'action à effectuer à la colonne j.

    start = emergence_simulation(P, blocks1, n_blocks, c_act, r_act);

    
    for (i = 0; i < n_blocks; i++) {    
      printf("%d ; %d \n", c_act[i], r_act[i]);
    }
    
    printf("nb d'actions déclanchées %d \n", start);

    // libération de la mémoire, fin du programme.
    free(blocks1);
    free(Q);
    free(where);
    free(c_act);
    free(r_act);
    //free(not_empty);
    spasm_csr_free(B);
    spasm_csr_free(A);
    uptri_free(P);
    exit(0);

    /* déterminer le nombre de diagonale non vide.
    n_diags = 0;
    for (i = 0; i < n_blocks; i++) {
      if(D[i] != -1) n_diags++;
    }
    */
    //printf("%d\n", n_diags);

    
    //for (k = 0; k < n_blocks; k++) {
    //printf("%d ; %d\n", k, last_b[k]+1);
    // }
     


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
