#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/*
 * Dans toute la suite, X[i:j] désigne les élements de X d'indice i <= ... < j,
 * C'est-à-dire de i inclus à j exclus. Il s'ensuit que i:j désigne l'intervalle [i; j[
 *
 * Cette notation bien pratique est héritée du langage Python.
 */



/****************** Structures ***********************/


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
  block_t *blocks; // coordonnées du du bloc corespondant(optionnel)
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
 * Type de donnée qui décrit le graphe "sympa" pour prévoir l'apparition des blocs.
 */
typedef struct {
  int entry; 
  int next; // entrée suivante sur la ligne.
  int col; // indice de colonne sur laquelle se trouve l'arrête
  //int diag; // indice de la diagonale correspondante
} edge_t;


/*
 * Liste chaînée d'entier.
 */
typedef struct list
{
  int val;
  struct list *next;
} list_t;


/*
 * Arbres binaire de recherche.
 */
typedef struct tree
{
  int val;
  struct tree *left;
  struct tree *right;
} tree_t;



/************************ allocation et libération de mémoire *************************/ 


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
uptri_t * uptri_alloc(int nzmax, int n, int with_blocks) {
  uptri_t *T;

  T = spasm_malloc(sizeof(uptri_t));

  T->nzmax = nzmax;
  T->n = n;
  T->d = spasm_malloc( (n+1) * sizeof(int));
  T->i = spasm_malloc(nzmax * sizeof(int));
  T->blocks = (with_blocks ? spasm_malloc(nzmax * sizeof(block_t)) : NULL);


  return T;
}

/*
 * libère la mémoire d'un uptri_t;
 */
void uptri_free(uptri_t *T) {
  if (T == NULL) return;

  free(T->d);
  free(T->i);
  free(T->blocks);
  free(T);
}


/*
 * vide une liste chaînée, libère la mémoire.
 */
void clear_list(list_t **list) {
  list_t *tmp;

  if(*list == NULL) return;

  tmp = *list;
  if(tmp->next) clear_list(&tmp->next);
  free(tmp);
  list = NULL;
}


/*
 * Vide un arbre binaire, libère la mémoire.
 */
void clear_tree(tree_t **tree) {
  tree_t *tmp;

  if(*tree == NULL) return;

  tmp = *tree;
  if(tmp->left != NULL) clear_tree(&tmp->left);
  if(tmp->right != NULL) clear_tree(&tmp->right);

  free(tmp);
  *tree = NULL;

}

/************************ fonctions diverses *******************************/

/* * * * * * * * * * * * * sur les blocs * * * * * * * * * * * */

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


/* * * * * * * * * * * * * sur les listes chaînées * * * * * * * * * * */

/*
 * Recherche l'élément a dans une liste chaînée d'entiers.
 * Retourne 1 si l'élément est dans la liste et 0 sinon.
 */
int exist_element(list_t *list, int a) {
  list_t *tmp;

  tmp = list;
  while (tmp) {
    if(tmp->val == a) return 1;
    tmp = tmp->next;
  }
  return 0;
}

/*
 * Ajouter un élément en tête de liste chaînée.
 */
list_t * new_element(list_t *list, int a) {
  list_t *new;
  new = spasm_malloc(sizeof(list_t));

  new->val = a;
  new->next = list;

  return new;
}

/* * * * * * * * * * * Autres * * * * * * * * * * */

/*
 * Recherche dichotomique dans une table tab, triée entre start et end.
 */
int dichotomie(int *tab, int start, int end, int val) {
  int found, mid;

  found = 0;

  while(!found && (end-start) > 1) {
    mid = (start + end)/2;
    found = (tab[mid] == val);
    if(tab[mid] > val) end = mid;
    else start = mid;
  }

  if (tab[start] == val) return start;

  else return -1;

}


/*********************** Décomposition par bloc *****************************/

/*
 * Stocke le L de la décomposition LU dans L.
 *  
 */
spasm_lu * sorted_submatrix_LU(const spasm *M, int a, int b, int c, int d, int *py) {
  spasm *C;
  int *p;
  spasm_lu *LU;

  //extrait la sous_matrice
  C =sorted_spasm_submatrix(M, a, c, b, d, py, SPASM_WITH_NUMERICAL_VALUES);
  if (spasm_nnz(C)==0) {
    spasm_csr_free(C);
   return NULL;
  }

  //calcule la décomposition LU
  p = spasm_cheap_pivots(C);
  LU = spasm_LU(C, p, SPASM_KEEP_L); // on garde L
  free(p);
 
  // libère la mémoire en trop.
  spasm_csr_free(C);

  return LU;
  
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
	    blocks[*start].r = 0; // initialisation.
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
  // spasm_lu **LU;

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

  // étape 2 : allouer la liste des blocs et du L correspondant

  blocks = spasm_malloc(sizeof(block_t) * k);
  *blocks_ptr = blocks;
  // LU = spasm_malloc(k * sizeof(spasm_lu *));
  // *LU_ptr = LU;
  
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

  // étape 5 : calculer les L et les rangs
  /* for (i = 0; i < k; i++) { */
  /*   LU[i] = sorted_submatrix_LU(M, blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1, py); */
  /*   blocks[i].r = LU[i]->U->n; */
   
  /* } */


  return k;
}

/*
 * Etant donné le rang des blocks diagonaux, Bi, d'une matrice M, triangulaire par
 * blocs, détermine pour chaque Bi le nombre de pivot qu'il reste à trouver sur
 * l'intervalle de ligne d'une part et l'intervalle de colonne d'autre part.
 */
void remaining_pivots_init(int *R, int *C, const block_t *blocks, int n_blocks ) {
  int i, a;
  for (i=0; i<n_blocks; i++) {

    // Intervalles des lignes :
    a = blocks[i].i0 + blocks[i].r; //<--- numéro de la ligne à partir de laquelle il n'y a plus de pivots
    R[i] = blocks[i].i1 - a; 

    // Intervalles des colonnes :
    a = blocks[i].j0 + blocks[i].r; //<--- numéro de la colonne à partir de laquelle il n'y a plus de pivots
    C[i] = blocks[i].j1 - a;
   
  }

}

/*
 * Etant donné un intervalle de ligne row et un bloc, met à jour le nombre de pivots
 * qu'il reste à trouver sur l'intervalle row.
 */
int row_pivot_update(int row, block_t block) {
  int a, b, val;

  a = block.i0 + block.r;
  b = block.i1 - a;
  val = (row > b) ? row - b : 0;

  return val;
}

/*
 * Etant donnée un intervalle de colonne col et un bloc, met à jour le nombre de pivots
 * qu'il reste à trouver sur l'intervalle col.
 */
int col_pivot_update(int col, block_t block) {
  int a, b, val;

  a = block.j0 + block.r;
  b = block.i1 - a;
  val = (col > b) ? col - b : 0;

  return val;
}

/*
 * Etant donnée une matrice diagonale par blocs, renvoie le
 * numéro de l'intervalle de colonne correspondant à la colonne j pour tout j.
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
 * Etant donné une matrice, le tableau de ses blocs diagonaux blocks,
 * la matrice de la position de ses blocs, un tableau Q
 * tel que Q[j] représente le numéro de l'intervalle de colonnes auquel 
 * j appartient, renvoie la liste des sous-matrices colonnes par 
 * colonnes (le premier élément de la liste est celui le plus près de la 
 * diagonale principale)
 */
void list_of_submatrices(spasm *A, spasm *B, block_t * blocks, int *Q, spasm_list **C){
  int Bn, k;
  int *Cm, *Cjstart;

  if(A == NULL || B == NULL){ // check inputs
    return;
  }
  Bn = B->n;
  assert(Bn = B->m); // B square.

 
  // get workspace :
  Cm = spasm_malloc(Bn * sizeof(int));
 Cjstart = spasm_malloc(Bn * sizeof(int));

  for(k = 0; k < Bn; k++){
    // if(blocks[k].j0 + blocks[k].r < blocks[k].j1){ // Il reste des pivots à trouver sur cet intervalle.
      Cm[k] = blocks[k].j1 - blocks[k].j0;
      Cjstart[k] = blocks[k].j0;
      assert(Cm[k] > 0);
      // }
  }

  for(k = 0; k < Bn; k++){
    spasm_list_of_submatrices_update(A, blocks[k].i0, blocks[k].i1, k, B, Q, Cm, Cjstart, C);

    }

  // free workspace :
  free(Cm);
  free(Cjstart);
  // free(w);

}

/******************** recherche d'intervalles "non complets" *******************/

/*
 * teste si un intervalle de ligne défini par le bloc diagonal "block"
 * est "complet" c-à-d si il y a un pivot sur chaque ligne de l'intervalle.
 *
 * renvoie 1 si c'est le cas et 0 sinon.
 */
int full_rows(block_t block, int r) {
  int t;

  t = block.i0 + r;
  assert(t <= block.i1);

  return (t == block.i1 ? 1 : 0);
}

/*
 * teste si un intervalle de colonne défini par un bloc diagonal "block" est 
 * "complet" c-à-d si il y a un pivot sur chaque ligne de l'intervalle.
 *
 * renvoie 1 si c'est le cas et 0 sinon.
 */
int full_columns(block_t block, int r) {
  int t;

  t = block.j0 + r;
  assert(t <= block.j1);

  return (t == block.j1 ? 1 : 0); 
}

/*
 * Donne la table des intervalles de lignes non complets.
 * la valeur renvoyée est le nombre d'entrée du tableau
 */
int rows_to_watch(block_t *blocks, int *r_tab, int n_blocks) {
  int count, i;

  count = 0;
    for(i = 0; i < n_blocks; i++) {
      if (!full_rows(blocks[i], blocks[i].r)) {
	r_tab[count] = i;
	count++;
      }
    }
  return count;
}

/********************* Recherche des blocs non vide ********************/


/*
 * Détermine le nombre de blocs non vide et les diagonales auquels ils appartiennent.
 *
 * Ici seul le dernier bloc observé à la diagonale k est stocké à D[k]. On n'a pas
 * d'information sur les autres blocs.
 *
 * Pour chaque intervalle de colonnes, détermine également le nombre de lignes non vides dans l'intervalle en question. 
 *
 * Utilise un tableau auxiliaire w, initialisé à -1 pour chaque intervalle de colonnes. si w[k] != -1 alors w[k] représente la dernière ligne rencontrée.
 *
 * La valeur renvoyée est le nombre de blocs non vides.
 * Le nombre de lignes non vides par intervalle de colonne est stocké dans le tableau n_rows.
 */
int count_non_empty_blocks_and_rows(const spasm *M, const block_t *blocks, int n_blocks, const int *Q, int *n_rows) {
  int i, l, p, j, *Mp, *Mj, n, c, k, fill, *D, *w;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  l = 0; // <--- numéro du bloc de ligne qu'on regarde.
  fill = 0;

  D = spasm_malloc(n_blocks * sizeof(int));
  w = spasm_malloc(n_blocks * sizeof(int));
  for (i = 0; i < n_blocks; i++) {
    D[i] = -1;
    w[i] = -1;
    n_rows[i] = 0;
  }

  for (i = 0; i < n; i++) {
    while (blocks[l].i1 <= i) l++;

    for (p = Mp[i]; p < Mp[i+1]; p++) {
      j = Mj[p];
      c = Q[j]; // <--- numéro de l'intervalle de colonne.
      if(w[c] != i){
	n_rows[c]++;
	w[c] = i;
      } 
      k = c - l; // <--- numéro de la diagonale à laquelle appartient l'entrée.
      if (D[k] != l) {
	fill++;
	D[k] = l; // <--- l remplace la valeur précédente de D[k].
      }
    }
  }
  free(D);
  free(w);
  return fill;
}

/*
 * Détermine la liste des emplacements des blocs non vides.
 * Pour chaque intervalle de colonnes, détermine la liste P des lignes non vides.
 * La valeur renvoyée est le nombre de blocs non vides.
 */
int non_empty_blocks_and_rows_list(const spasm *M, const block_t *blocks, int n_blocks, int * n_rows, const int *Q, blk_t *where, int ***P_pt) {
  int i, l, p, j, *Mp, *Mj, n, c, k, *D, count, **P, *w;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  l = 0; // <--- numéro du bloc de ligne qu'on regarde.
  //fill = 0;
  count = 0;

 P = spasm_malloc(n_blocks * sizeof(int*));
  *P_pt = P;

  //Get workspace.
  D = spasm_malloc(n_blocks * sizeof(int));
  w = spasm_malloc(n_blocks * sizeof(int));
  for (i = 0; i < n_blocks; i++) {
    D[i] = -1;
    P[i] = spasm_malloc(n_rows[i] * sizeof(int));
    n_rows[i] = 0; // réinitialiser le compteur à 0.
    w[i] = -1;
  }

  for (i = 0; i < n; i++) {
    while (blocks[l].i1 <= i) l++;

    for (p = Mp[i]; p < Mp[i+1]; p++) {
      j = Mj[p];
      c = Q[j]; //<--- numéro de l'intervalle de colonne.
      if(w[c] != i){
	P[c][n_rows[c]] = i;
	n_rows[c]++;
	w[c] = i;
      }
 
      k = c - l; // <--- numéro de la diagonale à laquelle appartient l'entrée.
      if (D[k] != l) {
	//fill++;
	D[k] = l; // <--- l remplace la valeur précédente de D[k].
	//if ((k == 0) || (blocks[c].r < blocks[c].j1 - blocks[c].j0)) {
	  where[count].c = c;
	  where[count].r = l;
	  count++;
	  // }
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
uptri_t * position_uptri(const blk_t *w, int n_blocks, int count, int with_blocks) {
  int k, *wd, *Bd, *Bi, *tmp, sum, p;
  uptri_t *B;

  /* Pour chaque bloc w[k] trouver sa diagonale. */

  wd = spasm_malloc(count * sizeof(int)); // Allocation mémoire diagonale

  for (k = 0; k < count; k++) {
    wd[k] = w[k].c - w[k].r; // <--- diagonale à laquelle appartient le bloc w[k]
  }

  /* allocation du résultat */
  B = uptri_alloc(count, n_blocks, with_blocks); // <--- On regarde les lignes

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




/************** Apparition de nouveaux blocs : prédiction de la structure *********************/

/*
 * Donne la transposée de la matrice de positions des blocs, à partir de la liste "where"
 * du nombre de blocs "count" et du nombres de blocs diagonaux "n_blocks".
 */
spasm * blocks_spasm(const blk_t *w, int n_blocks, int count, int prime, int transpose) {
  int k, *Tp, *Tj, *tmp, sum, p;
  spasm *T;

  /* allocation mémoire */
  T = spasm_csr_alloc(n_blocks, n_blocks, count, prime, 0);

  tmp = spasm_calloc(n_blocks, sizeof(int));

  Tp = T->p;
  Tj = T->j;

  /* compte les entrées sur les lignes de la matrice */
  for(k = 0; k < count; k++) {
    if(transpose) tmp[w[k].c]++;
    else tmp[w[k].r]++;
  }

  /* implémente les pointeurs sur les lignes de la matrice */
  sum = 0;
  for (k = 0; k < n_blocks; k++) {
    Tp[k] = sum;
    sum += tmp[k];
    tmp[k] = Tp[k];
  }
  Tp[n_blocks] = sum;

  /* trouve l'emplacement des entrées sur chaque ligne */
  for(k = 0; k < count; k++) {
    p = (transpose ? tmp[ w[k].c ]++ : tmp[w[k].r]++);
    Tj[p] = (transpose ? w[k].r : w[k].c);
  }

  /* libération de mémoire, retourner le résultat. */
  free(tmp);
  return T;
}



/*
 * On écrit la matrice d'adjacence du graphe "sympa".
 */
spasm * row_intersection_matrix(edge_t *row, int n_blocks, int count) {
  int sum, p, k, *Sp, *Sj, *Sx, *w;
  spasm *S;

  /* Allocate result */
  S = spasm_csr_alloc(n_blocks, n_blocks, count, 2147483647, 1);

  /*get workspace */
  w = spasm_calloc(n_blocks, sizeof(int));
  Sp = S->p;
  Sj = S->j;
  Sx = S->x;

  /* compute row counts */
  for (k = 0; k < count; k++) {
    w[ row[k].entry ]++;
  }

  /* compute row pointers */
  sum = 0;
  for (k = 0; k < n_blocks; k++) {
    Sp[k] = sum;
    sum += w[k];
    w[k] = Sp[k];
  }
  Sp[n_blocks] = sum;

  /* dispatch entries */
  for(k = 0; k < count; k++) {
    p = w[row[k].entry]++;
    Sj[p] = row[k].next;
    Sx[p] = row[k].col;
  }

  /* free workspace */
  free(w);
  return S;

}


/*
 * On fait pointer chaque entrée sur la suivante pour obtenir le graphe "sympa".
 */
spasm * row_intersection_graph(spasm *T, edge_t *rows, int n_blocks) {
  int k, i, *Tj, *Tp, n_entry, exist, j, jplus;
  list_t **pred; //liste de prédécesseurs du graphe.
  
  //get workspace
  pred = spasm_malloc(n_blocks * sizeof(list_t));

  for (i = 0; i < n_blocks; i++) {
    pred[i] = NULL;
  }

  Tp = T->p;
  Tj = T->j;

  n_entry = 0;

  for (i = 0; i < n_blocks; i++) {
   
    for (k = Tp[i]; k < Tp[i+1] - 1; k++) {
      j = Tj[k];
      jplus = Tj[k+1];

      exist = exist_element(pred[jplus], j);
      if (!exist) {
	rows[n_entry].entry = j;    // entrée courante.
	rows[n_entry].next = jplus; // entrée suivante sur la ligne.
	rows[n_entry].col = i;      // numéro de la colonne correspondante.
	n_entry ++;

	pred[jplus] = new_element(pred[jplus], j);
      }
    }
  }
  spasm *toto = row_intersection_matrix(rows, n_blocks, n_entry);


  // free workspace
  for (i = 0; i < n_blocks; i++) {
    clear_list(&pred[i]);
  }
  free(pred);

  return toto;
}

spasm * filled_structure(const spasm *A, const spasm *adjacency_graph) {
  spasm *B;
  int n, bnz, i, j, k, l, x_size, col, r, tmp_size;
  int *w, *x, *Bp, *Bj, *Ap, *Aj, *agj, *agp, *agx, *tmp_j, *tmp_p;

  // allouer résultat B (taille de A + taille de G)
  n = A->n;
  assert(adjacency_graph->n == n);

    /* get workspace */
  x = spasm_malloc(n * sizeof(int));
  w = spasm_malloc(n * sizeof(int));
  tmp_size = 4 * spasm_nnz(A);
  tmp_j = spasm_malloc(tmp_size * sizeof(int));
  tmp_p = spasm_malloc((n+1) * sizeof(int));

  /* allocate result */
  B = spasm_csr_alloc(n, n, 0, A->prime, 0);
  

  for (i = 0; i < n; i++) {
    /* clear workspace */
    x[i] = 0;
    w[i] = 0;
  }

  Ap = A->p;
  Aj = A->j;
  Bj = B->j;
  Bp = B->p;
  agj = adjacency_graph->j;
  agp = adjacency_graph->p;
  agx = adjacency_graph->x;


 
  // pour chaque ligne i en partant de la fin

  bnz = 0;
  for(i = n - 1; i >= 0; i--) {
    // invariant : pour toutes les lignes d'indice supérieur à i,
    //    * le début des indices de colonne est à tmp_p[i]
    //    * la taille de le k-ème ligne est dans Bp[k+1]

    //    dispatcher A[i] dans w et ajouter les entrées de A[i] dans x
    x_size = 0;
    for(k = Ap[i]; k < Ap[i+1]; k++) {
      assert(x_size < n);
      j = Aj[k];
      w[j] = 1;
      x[x_size] = j;
      x_size++;
    }


      //    pour chaque arête j --> i (arrive à la ligne i)
      // printf("%d\n", i);
      for(l = agp[i]; l < agp[i+1]; l++) {
	j = agj[l];
	col = agx[l];

	if (j <= i) {
	  printf("bug %d %d\n", i, j);
	  exit(1);
	  //      assert(j > i);
	}

	//        Dispatcher la ligne j de B à partir de la colonne c, dans w, et ajouter à x les nouvelles entrées
	for(k = tmp_p[j]; k < tmp_p[j] + Bp[j+1]; k++) {
	  r = tmp_j[k];
	  if (r >= col && w[r] == 0) {
	    w[r] = 1;
	    x[x_size] = r;
	    x_size++;
	  }
	}
      }

    //    si B n'est pas assez gros pour recevoir x, réallouer B 2x plus gros
    if (bnz + x_size > tmp_size) {
      //printf("tmp_size : %d ; %d\n", tmp_size, 2*tmp_size +n);
      tmp_size = 2 * tmp_size + n;
      tmp_j = realloc(tmp_j, tmp_size *(sizeof(int)));
      assert(tmp_j != NULL);
    }

    //    B[i] <-- x (copier les entrées ET ajuster pointeurs de ligne)
    //    Parcourir B[i] pour remettre w à zéro
    Bp[i+1] = x_size;
    tmp_p[i] = bnz;
    for(k = 0; k < x_size; k++) {
      j = x[k];
      tmp_j[bnz] = j;
      w[j] = 0;
      bnz++;
    }
  }

  //printf("nombre de blocs à alouer : %d\n", bnz);
  // printf("nombre de blocs initialement : %d\n", A->nzmax);
  assert(bnz >= A->nzmax);
  // Bp[i+1] = taille de la i-ème ligne
  // tmp_p[i] = position du début de la ième ligne dans tmp_j
  
  // somme-préfixe pour avoir les bons Bp
  Bp[0] = 0;
  for(i = 1; i < n+1; i++) {
    Bp[i] += Bp[i - 1];
  }

   
  // allouer Bj à la bonne taille et copier tmp_j dedans
  spasm_csr_realloc(B, bnz);
  Bj = B->j;

  for(i = 0; i < n; i++) {
    l = 0;
    for(k = Bp[i]; k < Bp[i+1]; k++) {
      Bj[k] = tmp_j[tmp_p[i] + l];
      l++;
    }
  }
  
  // libérer l'espace
  free(x);
  free(w);
  free(tmp_j);
  free(tmp_p);

  return B;
}


/*
 * réécrit la matrice de tous les blocs qui sont ou peuvent apparaître sous forme uptri_t
 * (vision diagonale par diagonale). Ne prend en compte que les blocs à partir de la première
 * diagonale supérieure et sur des intervalles de lignes et de colonnes "non complet".
 */
uptri_t * diagonal_structure(const spasm *B, int *r_tab, int n_rows) {
  int i, k, d, n_blocks, sum, px, py, *Bp, *Bj, *w, *Td, *Ti, nnz;
  uptri_t *T;

  n_blocks = B->n;
  Bp = B->p;
  Bj = B->j;

  assert(n_blocks >= n_rows);

  /* compute nnz, the number of entries in T*/

  nnz = 0;

  // Pour chaque intervalle de ligne "non complet" compter le nombre d'entrées moins celle sur la diagonale principale.
  for(k = 0; k < n_rows; k++) {
    i = r_tab[k];
    nnz += Bp[i+1] - Bp[i] - 1;  
  }

  /*Allocate result*/
  T = uptri_alloc(nnz, n_blocks, 0);
  Td = T->d;
  Ti = T->i;

  /* get workspace */
  w = spasm_calloc(n_blocks, sizeof(int));

  /* compute diagonal counts */
  for (k = 0; k < n_rows; k++) {
    i = r_tab[k];
    for(px = Bp[i]; px < Bp[i+1]; px++) {
      d = Bj[px] - i;
      w[d]++;
    }
  }        
  w[0] = 0; // on supprime les entrées sur la diagonale principale.

  /* compute diagonal pointers */
  sum = 0;
  for(d = 0; d < n_blocks; d++) {
    Td[d] = sum;
    sum += w[d];
    w[d] = Td[d];
  }
  Td[n_blocks] = sum;


  /* dispatch entries */
  for(k = 0; k < n_rows; k++) {
    i = r_tab[k];
    for(px = Bp[i]; px < Bp[i+1]; px++) {
      d = Bj[px] - i;
      // si on n'est pas sur la diagonale principale ajouter l'entrée
      if(d != 0) {
	py = w[d];
	Ti[py] = i;
	w[d]++;
      }
      //printf("Ti[%d] = %d \n", py, i);
    }
  }

  /* free memory */
  free(w);
  
  /* return */
  return T;

}

/****************** Traite les blocs sur les diag supérieures ********************/


/*
 * Extrait une sous-matrice M sur la première diagonale supérieure.
 * Calcule le produit des lignes de L^(-1) à partir de r (en tenant compte des
 * permutations).
 * retourne la sous-matrice de ce produit.
 */
spasm * lazy_submatrix_first_diag(const spasm *M, int *p, int r, spasm *L) {
  spasm *R;
  int l, l_new, i, px, *yi, ynz, n, Rm, nzmax, Rn, prime, *Rp, *Rj;
  spasm_GFp *y, *Rx;

  // vérifier les entrées.
  assert(r < M->n);
  assert(r < M->m);
  if(M == NULL || L == NULL) return NULL;
  

  //Allouer le résultat (nombre d'entrée dans entre i0 et i1).
  n = M->n;
  Rm = M->m;
  nzmax = M->nzmax;
  Rn = n - r;
  prime = M->prime;

  R = spasm_csr_alloc(Rn, Rm, nzmax, prime, SPASM_WITH_NUMERICAL_VALUES);
  Rp = R->p;
  Rj = R->j;
  Rx = R->x;

  //intialiser Rp.
  Rp[0] = 0;

  //Allouer mémoire de y et de yi et initialiser à 0.
  y = spasm_malloc(Rm * sizeof(spasm_GFp));
  yi = spasm_malloc(Rm * sizeof(int));
  spasm_vector_zero(y, Rm);
  spasm_vector_zero(yi, Rm);

  i = 0;

  // Pour k allant de r à n :
  for(l = r; l < n; l++) {
    //   effectuer les permutations de lignes.
    l_new = p[l];
    //   calculer le produit de la k-ième ligne de L^(-1) par M.
    //   le stocker dans y, et son support dans yi.
   
    ynz = spasm_inverse_and_product(L, M, l_new, y, yi, p);
    
    //   mettre à jour Rp.
    Rp[i+1] = Rp[i] + ynz;

    //   réallouer de la mémoire si besoin.
    if(nzmax < Rp[i+1]) {
      nzmax = 2 * nzmax + ynz;
      spasm_csr_realloc(R, nzmax);
      Rj = R->j;
      Rx = R->x;
    }
 
    //   ajouter les entrées dans Rj et Rx.
    for(px = 0; px < ynz; px++) {
      Rj[ Rp[i] + px ] = yi[px];
      Rx[ Rp[i] + px ] = y[ yi[px] ];
    }
    i++;
  }
  // Finaliser, ajuster le nombre d'entrée de R.
  spasm_csr_realloc(R, -1);
  
  // Libérer mémoire auxiliaire.
  free(y);
  free(yi);

  // Retourner R.
  return R;
}


/*
 * Extrait une sous-matrice, sur la première diagonale supérieure.
 * Calcule le produit des bonnes lignes de L^(-1)
 * retourne la sous-matrice de ce produit.
 */


/*
 * extrait une sous matrice sur une diagonale supérieure,
 * pour chaque ligne sur laquelle on a pas encore trouvé de pivot,
 * effectue la multiplication par le L correspondant, continue la décomposition LU
 * en utilisant cette ligne.
 *
 * la valeur renvoyée est le nombre de pivot supplémentaire trouvés.
 */
int upper_block_treatment(spasm_list *C, int k, int d, int ri, int **p, int *qinv, spasm_system **L, spasm *U) {
  spasm *R, *Lnew, *LM;
  int m, Rn, Un, r, top, i, deff, old_rank, unz, lnz, npiv, left;
  int *pnew, *xi, *ptmp, *Lp, *Up;
  spasm_GFp *x;

  /* fait les calculs paresseux. */

  LM = L[k]->M;
  old_rank = LM->m;
  assert(d>0); // d est la diagonale courrante supérieure à la diagonale principale.

  if(d == 1){  //<---- d représente ici la diagonale courante.
    spasm *M;
    assert(C->row = k); // Il y a une sous-matrice non vide sur la première diag supérieure.
    M = C->M;
    R = lazy_submatrix_first_diag(M, p[k], old_rank, LM);
  }
  else {
    // Pour l'instant retourne 0.
    return 0;
  }

  if(spasm_nnz(R)==0){
    return 0; //check matrix
  }

  /* Continue LU */

  assert(R->m == U->m);

  m = R->m; // nombre de colonne
  Rn = R->n; // nombre de ligne de R.
  Un = U->n; //nombre de ligne de U concaténé avec R.
  r = spasm_min(Un + Rn, m); // nombre de ligne au maximum dans le nouveau U

  assert(r > Un);

  deff = 0;
  npiv = 0;

  // Allocation de la mémoire de x la nouvelle ligne de la décomposition LU.
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(x, m);
  spasm_vector_zero(xi, 3*m);
  
  // Reallocation mémoire de U :

  spasm_csr_resize(U, r, m);
  unz = U->nzmax;
  Up = U->p;
  spasm_csr_realloc(U, 2 * U->nzmax + m);
 
  // Allocation mémoire de L :
  Lnew = spasm_csr_alloc(Rn, r, 2 * U->nzmax + m, R->prime, 1);
  lnz = 0;

  Lp = Lnew->p;
  // Allocation de la permutation pnew des lignes de Lnew.
  pnew = spasm_malloc(Rn * sizeof(int));

  // Initialisation :
  for(i = 0; i < Rn; i++){
    pnew[i] = i;
  }

  for(i = U->n; i < r; i++){
    Up[i] = 0;
  }

    // On prend les lignes une a une :
  for (i = 0; i < Rn; i++){
   
    Lp[i] = lnz;
    Up[i + Un - deff] = unz;

   //  On continue LU avec la ligne courrante :
    //-1- Réallocation de mémoire :
    if(unz + m > U->nzmax){
      spasm_csr_realloc(U, 2 * U->nzmax + m);
    }  
    if(lnz + m > Lnew->nzmax){
      spasm_csr_realloc(Lnew, 2 * Lnew->nzmax + m);
    }
    // -2- On résout le système triangulaire.
    top = spasm_sparse_forward_solve(U, R, i, xi, x, qinv);
    
    // -3- On trouve le pivot et on dispatche les coeffs dans la (r+i)-ème ligne de U et dans la i-ème ligne de L.

    npiv += spasm_find_pivot(xi, x, top, U, Lnew, &unz, &lnz, i, &deff, qinv, pnew, Rn);

    // spasm_vector_zero(x, m);
    // spasm_vector_zero(xi, 3*m);
    
  }
  // -4- On finalise Lnew et U.
  
  Up[i - deff] = unz;
  Lp[Rn] = lnz;

 //resize and realloc :
  spasm_csr_resize(U, i - deff, m);
  spasm_csr_resize(Lnew, Rn, Rn - deff);  
  spasm_csr_realloc(U, -1);
  spasm_csr_realloc(Lnew, -1);

  /* Mise à jour des structures de données pour les calculs paresseux suivants */

  // On trouve quel sera le prochain système "gauche"
  left = spasm_next_left_system(L, k, d-1);
  left = left - k;

 // On met à jour la liste L.
  L[k] = spasm_system_update(L[k], Lnew, pnew, old_rank, left, d);
 
  // On met à jour p.
  ptmp = spasm_malloc(Rn * sizeof(int)); // Allocation de mémoire temporaire.
  for(i = 0; i < Rn; i++){
    ptmp[i] = p[k][ri + pnew[i]];
  }
  for(i = 0; i < Rn; i++){
    p[k][ri + i] = ptmp[i];
  }
  free(ptmp);

  // libérer la mémoire
  spasm_csr_free(R);
  free(x);
  free(xi);
  // renvoie npiv, le nombre de nouveau pivots trouvé lors de cette décomposition. 
  return npiv;
} 


/*
 * Parcours la d-ième diagonale de la matrice par bloc, k >0.
 * Pour chaque bloc "intéressant" extrait la sous-matrice correspondante,
 * et effectue les oppérations du à la méthode paresseuse et calculer le rang.
 *
 * (Dans un premier temps, on ne regarde que la première diagonale supérieure.
 * On se contente ici de faire le produit par le L^(-1) du bloc diagonal à gauche du bloc
 * qu'on regarde actuellement et l'éliminition du bloc en même temps que la décomposition LU
 * La valeur renvoyée est le nombre de blocs sur la diagonale)
 */
int upper_research(spasm_list **C, const uptri_t *B, int d, const block_t *blocks, spasm_system **L, spasm **U, int **p, int **qinv, int *ri, int **npiv_ptr) {
  int px, k, i0, i1, j0, j1, rj, *Bd, *Bi, nbl, *n_piv;

  // vérifier les entrées.
  if(C == NULL || B == NULL) return -1;
  assert(d > 0 && d < B->n);

  Bd = B->d;
  Bi = B->i;
 
  // Compter le nombre d'élément dans la liste des blocs
  nbl = 0;
  for(px = Bd[d]; px < Bd[d+1]; px++) {
    k = Bi[px];
    i0 = blocks[k].i0;
    i1 = blocks[k].i1;
    j0 = blocks[k+d].j0;
    j1 = blocks[k+d].j1;
    rj = L[k]->M->m;
    
    if(i0 + ri[k] < i1 && j0 + rj < j1) {
      nbl++;
    }
  }

  //Allouer la mémoire de n_piv et de L_new.
  n_piv = spasm_malloc(nbl * sizeof(int));
  *npiv_ptr = n_piv;

  // parcourir la deuxième diagonale de la matrice des blocs.
  // Pour chaque bloc rencontré :
  nbl = 0;
  for(px = Bd[d]; px < Bd[d+1]; px++) {
    //   retrouver ses intervalles [i0 : i1] et [j0 : j1].
    k = Bi[px];
    i0 = blocks[k].i0;
    j0 = blocks[k+d].j0;
    i1 = blocks[k].i1;
    j1 = blocks[k+d].j1;
    rj = L[k]->M->m;

    //   trouver le rang du bloc diagonal en dessous. vérifier qu'il n'est pas complet en colonne.
    //   trouver le rang du bloc diagonal à sa gauche. vérifier qu'il n'est pas complet en ligne.
    if (i0 + ri[k] < i1 && j0 + rj < j1) {

    //   extraire la sous-matrice et calculer le produit.
      n_piv[nbl] = upper_block_treatment(C[k+d], k, d, ri[k], p, qinv[k+d], L, U[k+d]);

      //update ri[k]
      ri[k] += n_piv[nbl];

    nbl++;
    }
      
  }

  // retourner le nombre de bloc sur la diag.
  return nbl;

}


/**************** Fonction main *********************/

int main() {
  spasm_triplet *T;
  spasm *A, *B, *BP;
  spasm **U;  
  spasm_list **C;
  spasm_system **L;
  spasm_dm *x;
  spasm_lu **LU;
  int n_blocks, i, ne, nemax, prime, nbl;
  int *qinv, *Q, *ri, *n_rows; 
  int **p, *n_piv, **Uqinv, **rows_tab;
  block_t *blocks;
  blk_t *where;
  // spasm *Tr, *G;
  // int count, n_rows, Rank;
  // edge_t *rows;

  /* charge la matrice depuis l'entrée standard */

  T = spasm_load_sms(stdin, 42013);
  mem_alloc = 0;  
A = spasm_compress(T);

  spasm_triplet_free(T);

  prime = A->prime;

  /* met la matrice sous forme triangulaire par blocs
   * Calcule le rang des blocks de la première diagonale
   * garde en mémoire le L de la décomposition LU des blocs diagonaux */

  //calucle la décomposition de A
  x = spasm_dulmage_mendelsohn(A);

  // B = A permutée sous forme triangulaire par blocs
 
  qinv = spasm_pinv(x->DM->q, A->m);
  B = spasm_permute(A, x->DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  free(qinv);
  
  spasm_csr_free(A);
 

  /*
   * trouve les blocs qui composent B
   */

  n_blocks = block_list(B, x, &blocks);

  /*
   * Trouve les blocs non vides dans la matrice de départ.
   */

  Q = spasm_malloc(B->m * sizeof(int)); // Table qui à une colonne associe le numéro de son intervalle.
  column_diag_number(B, blocks, Q);

  n_rows = spasm_malloc(n_blocks * sizeof(int)); // Nombre de lignes non vide par intervalle de colonnes.

  nemax = count_non_empty_blocks_and_rows(B, blocks, n_blocks, Q, n_rows); // compte le nombre de blocs non vide.

  where = spasm_malloc(nemax * sizeof(blk_t));
  non_empty_blocks_and_rows_list(B, blocks, n_blocks, n_rows, Q, where, &rows_tab); //Donne la liste des blocs non vides.
  BP = blocks_spasm(where, n_blocks, nemax, prime, 0); //Matrice de la position des blocks.

 
  /*
   * Liste chainée de sous matrice par intervalle de colonnes.
   */
 /* mem_alloc = 0;  */
 /* C = spasm_malloc(n_blocks * sizeof(spasm_list*)); */
 /*  for(i = 0; i < n_blocks; i++){ */
 /*    C[i] = NULL; */
 /*  } */

 
 /*  list_of_submatrices(B, BP, blocks, Q, C); */
 
 /* LU = spasm_malloc(n_blocks * sizeof(spasm_lu*)); */

 /*  /\* */
 /*   * Traitement de la diagonale principale. */
 /*   *\/ */


 /*  for(i = 0; i < n_blocks; i++){ */
 /*    assert(C[i] != NULL); // Au moins une sous-matrice sur la diagonale. */
 /*    assert(C[i]->row == i); // le dernier bloc de l'intervalle de colonne est le bloc diagonal. */
   
 /*    assert(C[i]->M != NULL); */
 /*    assert(spasm_nnz(C[i]->M) != 0); */

  
 /*    LU[i] = spasm_LU(C[i]->M, SPASM_IDENTITY_PERMUTATION, SPASM_KEEP_L); */

 /*    C[i] = spasm_list_delete_first_matrix(C[i]); //<-- plus besoin du bloc diagonal. */
 /*    // C[i] = C[i]->up; //<-- le premier bloc de la liste est le premier au dessus de la diagonale. */

 /*    blocks[i].r = LU[i]->U->n; */
 
 /*    if(blocks[i].r == blocks[i].j1 - blocks[i].j0){ */
 /*      C[i] = spasm_list_free(C[i]); // plus de pivots à trouver sur l'intervalle de colonnes. */
 /*    } */

 /*  } */

 /*  /\* Affichage des résultats première diag *\/ */

 /*  for(i = 0; i < n_blocks; i++) { */
 /*    printf("%d : (%d, %d) -- (%d, %d), rank %d\n", i, blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1, blocks[i].r); */

 /*  } */

 /*  printf("--------------------------\n"); */

 /*  printf("blocs diagonaux : %d\n", n_blocks); */

 /*  printf("--------------------------\n"); */

 /*  /\* */
 /*   * Ecrit la matrice des positions des blocs sous forme "uptri_t" */
 /*   * (diagonale par diagonale) */
 /*   *\/ */

 /*  // Suprime les blocs situé sur un intervalle de colonnes non intéressant. */
  
 /*  ne = filled_blocks_list(B, blocks, n_blocks, Q, where); */

 /*  spasm_csr_free(B); */

 /*  // mise de la structure sous forme uptri_t pour l'avoir "diagonale par diagonale" */
 /*  uptri_t *DS = position_uptri(where, n_blocks, ne, 0); */
    
 /*  /\* */
 /*   * Initialisation des structures de données pour calculs paresseux. */
 /*   *\/   */

 /*  // récupérer le L et le U de la décomposition LU des blocs diagonaux. */
 /*  // Initialiser tableau ri et rj */
 /*  ri = spasm_malloc(n_blocks * sizeof(int)); */

 /*  L = spasm_malloc(n_blocks * sizeof(spasm_system)); */
 /*  U = spasm_malloc(n_blocks * sizeof(spasm*)); */
 /*  p = spasm_malloc(n_blocks * sizeof(int *)); */
 /*  Uqinv = spasm_malloc(n_blocks * sizeof(int *)); */
 /*  for(i = 0; i < n_blocks; i++) { */
 /*    L[i] = NULL; // Initialisation de la liste chainée. */
 /*    U[i] = LU[i]->U; */
 /*    p[i] = LU[i]->p; */
 /*    Uqinv[i] = LU[i]->qinv; */
 /*    L[i] = spasm_system_update(L[i], LU[i]->L, LU[i]->p, 0, 0, 0); */
 /*    ri[i] = U[i]->n; */
 /*  } */

 /*  /\* */
 /*   * Traitement de la première diagonale supérieure. */
 /*   * Pas encore de remplissage. */
 /*   *\/ */

 /*  nbl = upper_research(C, DS, 1, blocks, L, U, p, Uqinv, ri, &n_piv); */

 /*  printf("nombre de blocs à traiter sur la première diagonale supérieure : %d\n", nbl); */

 /*  for(i = 0; i < nbl; i++) { */
 /*    printf("1ere diag sup : block %d : rank : %d \n", i, n_piv[i]); */
 /*   } */

 /*  printf("---------------------\n"); */

 
/*   // remplissage des autres diagonales. */
/*   Tr = blocks_spasm(where, n_blocks, count, prime, 1); */
/*   spasm_row_entries_sort(Tr, 0); */

/*   rows = spasm_malloc((spasm_nnz(Tr) - n_blocks) * sizeof(edge_t)); */

/*   spasm *row_inter = row_intersection_graph(Tr, rows, n_blocks); */
/*   spasm * blocks_mat = blocks_spasm(where, n_blocks, count, prime, 0); */
/*   uptri_t *FS; */
 
/*   G = filled_structure(blocks_mat, row_inter); */

/*   int * r_tab = spasm_malloc(n_blocks *sizeof(int)); */

/*   n_rows = rows_to_watch(blocks0, r_tab, n_blocks); */
/*   FS = diagonal_structure(G, r_tab, n_rows); */

/*   free(r_tab); */
/*   free(rows);  */
/*   spasm_csr_free(Tr); */
/*   spasm_csr_free(row_inter); */
/*   spasm_csr_free(blocks_mat); */
/*   spasm_csr_free(G); */

/*   printf("nombre de blocs intéressants au total : %d\n", FS->nzmax + n_blocks); */
/*   printf("---------------------\n"); */

/*   // essai sur la deuxième diag sup : */
/*   // free(n_piv); */

/*   // nbl = upper_research(M, DS, 2, blocks0, L, U, LUp, LUqinv, ri, &n_piv);  */

/*   /\* printf("nombre de blocs sur 2 diag : %d\n", nbl); *\/ */

/*   /\* for(i = 0; i < nbl; i++){ *\/ */
/*   /\*   printf("i : %d nombre de pivots : %d\n", i, n_piv[i]); *\/ */
/*   /\* } *\/ */

  // libération de la mémoire, fin du programme.
  for(i = 0; i < n_blocks; i++) {
    /* L[i] = spasm_system_clear(L[i]); */
    /* spasm_csr_free(U[i]); */
    /* free(Uqinv[i]); */
    /* free(LU[i]); */
    /* spasm_list_free(C[i]); */
    free(rows_tab[i]);
 }
  
 
  free(blocks);
  free(n_rows);
  free(rows_tab);
  /* free(ri); */
  /* free(LU); */
  /* free(L); */
  /* free(U); */
  /* free(p); */
  /* free(Uqinv); */
  free(Q);
  free(where);
  // spasm_csr_free(BP);
  //free(n_piv);
  // uptri_free(DS);
  // uptri_free(FS);
 
  return 0;
}
