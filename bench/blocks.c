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
 * Type de donnée qui décrit le graphe "sympa" pour prévoir l'apparition des blocs.
 */
typedef struct {
  int entry; 
  int next; // entrée suivante sur la ligne.
  int col; // indice de colonne sur laquelle se trouve l'arrête
  //int diag; // indice de la diagonale correspondante
} edge_t;


/*
 * Type de donnée qui définie les actions prévues en une ligne/colonne.
 * Utilise les listes simplement chaînées
 *
typedef struct action 
{
  int act; // élément de la liste chaînée
  struct action *next; // pointeur sur l'élément suivant
} action_t; 

*/

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
 * Renvoie le rang de M[a:c, b:d]. 
 */
spasm_lu * submatrix_LU(const spasm *M, int a, int b, int c, int d) {
  spasm *C;
  //int *p;
  spasm_lu *LU;

  // extrait la sous-matrice
  C = spasm_submatrix(M, a, c, b, d, SPASM_WITH_NUMERICAL_VALUES);
  if (spasm_nnz(C) == 0) {
    spasm_csr_free(C);
    return NULL;
  }

  // calcule la décomposition LU. 
  // p = spasm_cheap_pivots(C); 
  LU = spasm_LU(C, SPASM_IDENTITY_PERMUTATION, SPASM_KEEP_L); // on garde L
  // free(p);


  // libère la sous-matrice et la mémoire dont on n'a plus besoin.
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
int block_list(const spasm *M, const spasm_dm *DM, block_t **blocks_ptr, spasm_lu ***LU_ptr) {
  int i, k;
  block_t *blocks;
  spasm_lu **LU;
  
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
  LU = spasm_malloc(k * sizeof(spasm_lu *));
  *LU_ptr = LU;
  
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
    LU[i] = submatrix_LU(M, blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1);
    blocks[i].r = LU[i]->U->n;
  }

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
 * numéro du bloc diagonal correspondant à la colonne j pour tout j.
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
uptri_t * position_uptri(const blk_t *w, int n_blocks, int count) {
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
uptri_t * final_structure_uptri(const spasm *B, int *r_tab, int n_rows) {
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
  T = uptri_alloc(nnz, n_blocks, 1, 0);
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

/* /\* */
/*  * met à jour les tables des pivots à chercher pour la diagonale diag. */
/*  *\/ */
/* void update_remaining_pivots(int diag, const spasm *M, const uptri_t *B, int *R, int *C, const block_t *blocks) { */
/*   int k, i, j, *Bd, *Bi, nbl, l; */
/*   block_t *current; */

/*   Bd = B->d; */
/*   Bi = B->i; */

/*   nbl = Bd[diag+1] - Bd[diag]; // <--- nombre d'entrée sur la diagonale k. */

/*   current = spasm_calloc(nbl, sizeof(block_t)); */
/*   l = 0; */

/*   // Parcours des entrées de la diagonale diag : */
/*   for(k = Bd[diag]; k < Bd[diag +1]; k++) { */
/*     i = Bi[k]; */
/*     j = i + diag; */
   
/*     //remplissage des données des blocks rencontrés sur la diagonales. */
/*     if(R[i] > 0 && C[j] > 0) { */
/*       current[l].i0 = blocks[i].i0; */
/*       current[l].i1 = blocks[i].i1; */
/*       current[l].j0 = blocks[j].j0; */
/*       current[l].j1 = blocks[j].j1; */
/*       current[l].r = submatrix_rank(M, current[l].i0, current[l].j0, current[l].i1, current[l].j1); */
    
/*       //mise à jour des tables R et C du nombre de pivots restants en ligne et en colonne. */
/*       R[i] = row_pivot_update(R[i], current[l]); */
/*       C[j] = col_pivot_update(C[j], current[l]); */

/*     } */
/*     l++; */
/*   } */
/*   free(current); */
/* } */

/* /\* */
/*  * Parcourt la matrice diagonale par diagonale et compte le nombre de pivot */
/*  * qu'il faut encore trouver par intervalle de lignes et de colonnes sur les diagonales */
/*  * supérieures. */
/*  * Le programme s'arrête quand on a trouvé tous les pivots. */
/*  * La valeur renvoyée est le numéro de la diagonale où le programme s'est arrêté. */
/*  *\/ */
/* int last_diag_estimation(const spasm *M, const uptri_t *B, const block_t *blocks, int n_blocks) { */
/*   int k, diag, *R, *C, Mn, Mm, nbl, start; */
/*   char matrix_type; */

/*   Mn = M->n; */
/*   Mm = M->m; */

/*   diag = 0; */

/*   matrix_type = (Mn > Mm) ? 'V' : 'H'; // regarde si la matrice est verticale ou horizontale. */

/*   //Allocation mémoire des tables R et C du nombre de pivots qu'il reste à trouver. */
/*   R = spasm_malloc(n_blocks * sizeof(int)); */
/*   C = spasm_malloc(n_blocks * sizeof(int)); */

/*   // Initialisations des tables R et C. */
/*   remaining_pivots_init(R, C, blocks, n_blocks); */

/*   switch (matrix_type) { */
/*   case 'H' : */
/*     k = 0;// indicateur qui parcours les intervalles de lignes */
/*     nbl = n_blocks- 1; */

/*     while(k < nbl) { */
/*       for ( ; k < nbl && R[k] == 0; k++); // parcourt R jusqu'à ce qu'on trouve une entrée non nulle. */
/*       diag++; */
/*       update_remaining_pivots(diag, M, B, R, C, blocks); */
/*       nbl--; */
/*     } */
/*     break; */
/*   case 'V' : */
/*     k = n_blocks - 1; */
/*     start = 0; */

/*     while(k > start) { */
/*       for( ; k > start && C[k] == 0; k--); // parcourt C à l'enver jusqu'à ce qu'on trouve une entrée non nulle. */
/*       diag++; */
/*       update_remaining_pivots(diag, M, B, R, C, blocks); */
/*       start++; */
/*     } */
/*     break; */
/*   default : */
/*     printf("Matrix error. \n"); */
/*     diag = -1; */
/*     break; */
/*   } */

/*   free(R); */
/*   free(C); */
/*   return diag; */
  
/* } */


/**************** Fonction main *********************/

int main() {
  spasm_triplet *T;
  spasm *A, *B;
  spasm_dm *x;
  spasm_lu **LU;
  int n_blocks, i, *qinv, rank;
  block_t *blocks;
  spasm *Tr, *G;
  int count, n_rows, *r_tab;
  edge_t *rows;

  //interval_t *R, *C;

  // charge la matrice depuis l'entrée standard
  T = spasm_load_sms(stdin, 42013);
  A = spasm_compress(T);
  spasm_triplet_free(T);

  //#ifdef BAAAAD

  // calcule la décomposition de A
  x = spasm_dulmage_mendelsohn(A);

  // B = A permutée sous forme triangulaire par blocs
  qinv = spasm_pinv(x->DM->q, A->m);
  B = spasm_permute(A, x->DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  free(qinv);
  //spasm_save_csr(stdout, B);

 
  // calcule la liste des blocs
  n_blocks = block_list(B, x, &blocks, &LU);
  rank = 0;
  for(i = 0; i < n_blocks; i++) {
    printf("%d : (%d, %d) -- (%d, %d), rank %d\n", i, blocks[i].i0, blocks[i].j0, blocks[i].i1, blocks[i].j1, blocks[i].r);
    rank += blocks[i].r;
  }
  printf("blocs diagonaux : %d\n", n_blocks);
  printf("borne inf sur le rang : %d\n", rank);

  free(x);

  int *Q, fill, start, last_diag, entries;
  blk_t *where;
   

  // allocation mémoire de not_empty et Q
  //not_empty = malloc(n_blocks * sizeof(int));
   Q = malloc(B->m * sizeof(int));

  // trouver le numéro de l'intervalle auquel appartient une colonne.
  column_diag_number(B, blocks, Q);
  printf("-------------------------\n");

  // trouver le nombre de diagonales ayant au moins un bloc non vide.
  // n_diags = diag_count(B, blocks1, Q, not_empty);

  //printf("%d\n", n_diags);


  fill = count_filled_blocks(B, blocks, n_blocks, Q); // <--- nombre total de blocs non vide.

  printf("nombre total de blocks non-vide : %d\n", fill);

  //allocation de mémoire de where.
  where = malloc(fill * sizeof(blk_t));
  for (i = 0; i < fill; i++) {
    where[i].c = -1;
    where[i].r = -1;
  }
    
  count = filled_blocks_list(B, blocks, n_blocks, Q, where);
  printf("nombre de blocs structure initiale %d\n", count);

  Tr = blocks_spasm(where, n_blocks, count, B->prime, 1);
  //#else
  // for debugging purposes
  //Tr = spasm_transpose(A, 0);
  //#endif

 printf("----------------------\n");

 /* remplissage : on oublie dans un premier temps */

  /* spasm_row_entries_sort(Tr, 0); */
  /* rows = spasm_malloc((spasm_nnz(Tr) - n_blocks) * sizeof(edge_t)); */

  
  /* spasm *row_inter = row_intersection_graph(Tr, rows, n_blocks); */
  /* spasm * blocks_mat = blocks_spasm(where, n_blocks, count, B->prime, 0); */
  /* uptri_t *FS; */

  

  /* printf("nombre d'arêtes graphe d'adjacence : %d\n", spasm_nnz(row_inter)); */
  /*  G = filled_structure(blocks_mat, row_inter); */
  
  /* printf("nombre de blocs dus au remplissage : %d\n", G->nzmax); */

  /* r_tab = spasm_malloc(n_blocks *sizeof(int));   */

  /* n_rows = rows_to_watch(blocks, r_tab, n_blocks); */
  /* FS = final_structure_uptri(G, r_tab, n_rows); */

  /* printf("nombre de blocs intéressants au total : %d\n", FS->nzmax + n_blocks); */

  /* libération de la mémoire, fin du programme. */

  for(i = 0; i < n_blocks; i++){
    spasm_free_LU(LU[i]);
  }

  free(Q);
  free(LU);
  //free(rows);
  //free(r_tab);
  free(blocks);
  spasm_csr_free(B);
  spasm_csr_free(A);
  spasm_csr_free(Tr);
  // spasm_csr_free(G);
  //spasm_csr_free(row_inter);
  // spasm_csr_free(blocks_mat);
  // uptri_free(FS);
  return 0;
}
