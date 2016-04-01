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
  int r;  // le nombre de pivot sur l'intervalle [i0; i1[
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
	  blocks[*start].r = 0; // initialement aucun pivots trouvé sur l'interval de ligne *start.
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
    while (blocks[l].i1 <= i) {
      l++;
    }

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
 * Détermine la liste des emplacements des blocs non vides.
 */
blk_t * filled_blocks_list(const spasm *M, const block_t *blocks, int n_blocks, const int *Q, int fill) {
  int i, l, p, j, *Mp, *Mj, n, c, k, *D, count;
  blk_t * where;

  Mp = M->p;
  Mj = M->j;
  n = M->n;
  l = 0; // <--- numéro du bloc de ligne qu'on regarde.
  count = 0; // <--- nombre de blocs qu'on a ajouté à la liste

  // Initialisation de la liste de blocs. 
  where = spasm_malloc(fill * sizeof(blk_t));
  for(i = 0; i < n_blocks; i++) {
    where[i].c = -1;
    where[i].r = -1;
  }

  // D[i] =  numéro du dernier intervalle de ligne rencontré sur la diagonale D
  D = spasm_malloc(n_blocks * sizeof(int));
  for (i = 0; i < n_blocks; i++) {
    D[i] = -1;
  }

  for (i = 0; i < n; i++) {
    while (blocks[l].i1 <= i) l++; // l = numéro de l'intervalle de lignes auquel i appartient. 

    for (p = Mp[i]; p < Mp[i+1]; p++) {
      j = Mj[p];
      c = Q[j]; 
      k = c - l; // <--- numéro de la diagonale à laquelle appartient l'entrée.
      if (D[k] != l) { 
	// j est la première entrée observée sur le bloc (l, c).
        D[k] = l; // <--- l remplace la valeur précédente de D[k].
	where[count].c = c;
	where[count].r = l;
	count++;
      }
    }
  }
  free(D);
  return where;
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
 * CHARLES THINKS : est-ce qu'on ne pourrait pas simplement utiliser la fonction qui convertit une matrice triplet en CSR ?
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
 * (vision diagonale par diagonale).
 */
uptri_t * final_structure_uptri(const spasm *B) {
  int k, d, n_blocks, sum, px, py, *Bp, *Bj, *w, *Td, *Ti, nnz;
  uptri_t *T;

  n_blocks = B->n;
  Bp = B->p;
  Bj = B->j;

  /* compute nnz, the number of entries in T*/

  nnz = 0;

  // Pour chaque intervalle de ligne "non complet" compter le nombre d'entrées.
  for(k = 0; k < n_blocks; k++) {
    nnz += Bp[k+1] - Bp[k];  
  }

  // printf("nnz : %d\n", nnz);
  /*Allocate result*/
  T = uptri_alloc(nnz, n_blocks, 1, 0);
  Td = T->d;
  Ti = T->i;

  /* get workspace */
  w = spasm_calloc(n_blocks, sizeof(int));

  /* compute diagonal counts */
  for (k = 0; k < n_blocks; k++) {
    for(px = Bp[k]; px < Bp[k+1]; px++) {
      d = Bj[px] - k;
      w[d]++;
    }
  }        

  /* compute diagonal pointers */
  sum = 0;
  for(d = 0; d < n_blocks; d++) {
    Td[d] = sum;
    sum += w[d];
    w[d] = Td[d];
  }
  Td[n_blocks] = sum;

  /* dispatch entries */
  for(k = 0; k < n_blocks; k++) {
    for(px = Bp[k]; px < Bp[k+1]; px++) {
      d = Bj[px] - k;
      // si on n'est pas sur la diagonale principale ajouter l'entrée
      //if(d != 0) {
  py = w[d];
  Ti[py] = k;
  w[d]++;
  // }
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
 * fonction qui prend en entrée un intervalle de n lignes, 
 * un nombre de pivot déjà trouvé sur
 * cet interval n_piv, une permutation p qui met les lignes avec
 * les pivots en "haut" de l'intervalle.
 *
 * Pour chaque ligne sur lesquelles il reste des pivots à trouver,
 * effectue les calculs paresseux, et "fait entrer" la ligne dans la
 * décomposition LU de la bonne matrice U.
 * Met à jour n_piv et la permutation p.
 *
 * Met à jour ln_ptr pointeurs sur le nombre 
 * de ligne actuel de L
 *
 * renvoie n_piv.
 */
int upper_block_research(super_spasm *CS, super_list *List, block_t block, super_spasm *L, super_spasm *U, int *p, int *qinv, int *ln_ptr, int *un_ptr, int *lnz, int *unz){
  int i, j, i_new, tmp_start, m, ynz, top, n, *xi, *yi, ln, un, found, deff, n_piv, *tmp, *Up, *Lp, re_ord;
  spasm_GFp *x, *y;

  //check inputs :
  assert(U != NULL);
  assert(L != NULL);
  assert(p != NULL);
  assert(qinv != NULL);

  m = CS->M->m;
  // printf("m : %d et %d\n", U->M->m, m);
  assert(U->M->m == m);
  Up = U->M->p;
  Lp = L->M->p;

  un = *un_ptr;
  ln = *ln_ptr;
  deff = 0;
  n_piv = 0;
  tmp_start = block.i0 + block.r;
  n = block.i1 - (tmp_start);

  /* get worspace */
  y = spasm_malloc(m * sizeof(spasm_GFp));
  yi = spasm_malloc(m * sizeof(int));
  x = spasm_malloc(m * sizeof(spasm_GFp));
  xi = spasm_malloc(3 * m * sizeof(int));
  tmp = spasm_malloc(n * sizeof(int));

  /* main loop : search pivots */
  for(i = block.i0 + block.r; i < block.i1; i++){
    i_new = p[i]; //<-- "good" row.

    //printf("i = %d\n", i);
    /*check permutation */
    if(i_new >= block.i1 || i_new < block.i0){
      printf("Permutation error \n");
      exit(0);
    }

    /* clear workspace */
    for(j = 0; j < m; j++){
      y[j] = 0;
      yi[j] = 0;
      x[j] = 0;
    }
    for(j = 0; j < 3*m; j++){
      xi[j] = 0;
    }
  
    Up[un] = (*unz);
    Lp[ln] = (*lnz);

    //Lazy computation :

    ynz = super_spasm_lazy(CS, List, i_new, y, yi); 

    /* continue LU */
    // not enough room : realloc
    if ((*lnz) + m > L->M->nzmax) {
      spasm_csr_realloc(L->M, 2 * L->M->nzmax + m);
    }
    if ((*unz) + m > U->M->nzmax) {
      spasm_csr_realloc(U->M, 2 * U->M->nzmax + m);
    }
   
    // solve system x * L = y :

    top = spasm_sparse_forward_solve_scat(U->M, y, yi, ynz, xi, x, qinv);
    // find pivot (if there is one) and dispatch x :
    found = super_spasm_find_pivot(xi, x, top, U, L, unz, lnz, ln, un, i_new, qinv);
    // update npiv, ln, un tmp:
    n_piv += found;
    ln++;
    un += found;
    if(found){
      tmp[i - tmp_start - deff] = i_new;
     
    }
    else{
      deff++;
      tmp[n - deff] = i_new;
    }
  }

  // vérifier qu'il n'y a pas d'incohérence.
  assert(tmp_start + n_piv <= block.i1);

  // mettre à jour la permutation p en la réordonnant :
  for(i = tmp_start; i < tmp_start + n_piv; i++){
    p[i] = tmp[i - tmp_start]; // mettre les nouveaux pivots à la suite des autres.
  }
  re_ord = (block.i1 - tmp_start) - 1;
  for(i = tmp_start + n_piv; i < block.i1; i++){
    p[i] = tmp[re_ord];
    re_ord--;
  }

  /*free workspace */
  free(tmp);
  free(x);
  free(y);
  free(xi);
  free(yi);

  /* update ln_ptr and un_ptr */
  *ln_ptr = ln;
  *un_ptr = un;

  return n_piv;
}


/**************** Fonction main *********************/

int main() {
  spasm_triplet *T;
  spasm *A, *B;
  super_spasm **CS, **U, **L;
  super_list *L_list = NULL; 
  spasm_dm *x;
  int n_blocks, i, k, *qinv, rank, prime, **Qinv, *P, n_big, m, nzmax, diag, pt;
  block_t *blocks;
  spasm *Tr, *G;
  edge_t *rows;

  /* --------- MISE EN PLACE DES STRUCTURES ----------*/

  /* Ecriture de la matrice par bloc */

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
  //spasm_save_csr(stdout, B);

  /* Vrai rang méthode LU pour vérifier les résultats */
  spasm_lu *LU = spasm_LU(A, NULL, 0);
  printf("vrai rang : %d\n", LU->U->n);
  spasm_csr_free(A);

  // calcule la liste des blocs diagonaux
  n_blocks = block_list(B, x, &blocks);

  free(x); // le vieil objet DM : exit

  int *J, fill;
  blk_t *where;
   
  /*---------------------------------------------------*/

  /*  Obtention du motif des blocs de la matrice B */

  // J[j] = numéro de l'intervalle auquel appartient la colonne j.
  J = malloc(B->m * sizeof(int));
  column_diag_number(B, blocks, J);

  fill = count_filled_blocks(B, blocks, n_blocks, J); // <--- nombre total de blocs non vide.

  // remplit where en donnant l'intervalle de ligne (where.r) et l'intervalle de colonne (where. c) de tous les blocs non vides.
  where = filled_blocks_list(B, blocks, n_blocks, J, fill); 

  // calcule la matrice des blocks
  Tr = blocks_spasm(where, n_blocks, fill, B->prime, 1);
  
  /*-------------------------------------------------*/

  /* Remplissage potentiel de la matrice de départ */

  spasm_row_entries_sort(Tr, 0);

  rows = spasm_malloc((spasm_nnz(Tr) - n_blocks) * sizeof(edge_t));
  
  spasm *row_inter = row_intersection_graph(Tr, rows, n_blocks);
  spasm * blocks_mat = blocks_spasm(where, n_blocks, fill, B->prime, 0);

  free(where);
  uptri_t *FS;

  G = filled_structure(blocks_mat, row_inter); // matrice "spasm" du remplissage
  

  FS = final_structure_uptri(G); // matrice "diag par diag" du remplissage.

  spasm_csr_free(G);
  spasm_csr_free(Tr);
  spasm_csr_free(row_inter);
  spasm_csr_free(blocks_mat);
 
  /* --------------------------------------------------- */


  /* Découpage de B en tranche de colonnes */
  int *first_col = spasm_malloc((n_blocks + 1) * sizeof(int));
  for(i = 0; i < n_blocks; i++){
    first_col[i] = blocks[i].j0;
  }
  first_col[n_blocks] = B->m;

  CS = super_spasm_column_slices(B, J, first_col, n_blocks, 1);
  n_big = B->n;
  m = B->m;
  nzmax = B->nzmax;
  prime = B->prime;

  free(first_col);
  spasm_csr_free(B);



  /* -------------------------------------------------------- */
  

  /* ----------- INTIALISATTION  -------------- */

  /* Get workspace */

  int *unz, *un, lnz, ln;
  U = spasm_malloc(n_blocks * sizeof(super_spasm*));
  L = spasm_malloc(n_blocks * sizeof(super_spasm*));
  L_list = NULL; //<-- initialement la liste des L est vide.
  Qinv = spasm_malloc(n_blocks * sizeof(int*));
  P = spasm_malloc(n_big * sizeof(int));
  unz = spasm_malloc(n_blocks * sizeof(int));
  un = spasm_malloc(n_blocks * sizeof(int));


  int U_nmax;
  for(k = 0; k < n_blocks; k++){
    U_nmax = spasm_min(CS[k]->M->m, n_big); // borne sup sur le nombre de pivots qu'on peut trouver dans l'intervalle de cols J_k.
    Qinv[k] = spasm_malloc(CS[k]->M->m * sizeof(int));
    U[k] = super_spasm_alloc(n_big, U_nmax, CS[k]->M->m, 2 * (U_nmax + CS[k]->M->m), prime, 1);
  }

  /* ------------------------------------------ */

  /* Initialize LU workspace */

  for(k = 0; k < n_blocks; k++){
    unz[k] = 0; // U initialement vide. 
    un[k] = 0; // On se place sur la ligne 0 de U->M

    for(i = 0; i < CS[k]->M->m; i++){
      Qinv[k][i] = -1; // pas encore de pivots trouvés.
    }
  }
  for(i = 0; i < n_big; i++){
    P[i] = i; // pas encore de pivots trouvés.
  }

  /* -------------------------------------------*/
  int piv_i = n_big; // nombre de lignes sur lesquels il reste (potentiellement) des pivot.
  diag = 0;
  int n_piv, last_block_rows, piv_on_diag = 0;
  int J_good = n_blocks; // nombre d'intervalle de colonnes sur lesquels il reste (potentiellement) des pivots. 

  /* MAIN LOOP : RECHERCHER LES PIVOTS EN REMONTANT LES DIAGONALES */
 
  while(piv_i > 0 && J_good >0){ //<--- tant qu'il reste des pivots potentiels
    /* allouer le L de la diagonale */
    L[diag] = super_spasm_alloc(n_big, piv_i, m, nzmax, prime, 1);
    /* initialisation */
    ln = 0; // On se place sur la ligne 0 de L
    lnz = 0; // L est pour l'instant vide.

    /* boucle sur tous les blocks (k, k+d) à regarder */
    for(pt = FS->d[diag]; pt < FS->d[diag + 1]; pt ++){
      k = FS->i[pt];

      // Si tous les pivots de l'intervalle de ligne I_k on déjà été trouvé on passe au suivant.
      if(blocks[k].i0 + blocks[k].r == blocks[k].i1){
      	continue;
      }

      //traitement du block k:
      n_piv = upper_block_research(CS[k+diag], L_list, blocks[k], L[diag], U[k+ diag], P, Qinv[k + diag], &ln, &un[k + diag], &lnz, &unz[k + diag]);


      blocks[k].r += n_piv;
      piv_i = piv_i - n_piv;
      piv_on_diag += n_piv;

      // Vérifier que le nombre de pivot trouvé est possible :
      assert(blocks[k].i0 + blocks[k].r <= blocks[k].i1);
      assert(un[k + diag] <= U[k+diag]->M->n);

      //Si U[k + diag] est terminé le finaliser et mettre à jour J_good:
      if(un[k + diag] == U[k + diag]->M->n){
	U[k + diag]->M->p[un[k + diag]] = unz[k + diag];
	spasm_csr_resize(U[k + diag]->M, un[k + diag], U[k + diag]->M->m);
	spasm_csr_realloc(U[k + diag]->M, -1);
	J_good--; //plus de pivots à chercher sur cet intervalle de colonnes.
      }

    }
    //Finaliser L[diag]:
    L[diag]->M->p[ln] = lnz;
    spasm_csr_resize(L[diag]->M, ln, n_big);
    spasm_csr_realloc(L[diag]->M, -1);

    //Si U[diag] pas encore finalisée, finaliser U[diag], mettre à jour J_good:
    if(un[diag] != U[diag]->M->n){
	U[diag]->M->p[un[diag]] = unz[diag];
	spasm_csr_resize(U[diag]->M, un[diag], U[diag]->M->m);
	spasm_csr_realloc(U[diag]->M, -1);
	J_good--; // on ne regarde plus l'intervalle de colonne le plus à gauche.
    }

    // mettre à jour la liste chainée L et diag :
    L_list = super_list_update(L_list, L[diag]);
    diag++;

    // mettre à jour piv_i
    last_block_rows = blocks[n_blocks - diag].i1 - (blocks[n_blocks - diag].i0 + blocks[n_blocks - diag].r);

    piv_i = piv_i - last_block_rows; // on ne recherche plus de pivot sur ces lignes.

  }


  /* -------- FINALISATION DU PROGRAMME -------------------*/

  /* printf("U :\n"); */
  /* for(k = 0; k < n_blocks; k++){ */
  /*   spasm_save_csr(stdout, U[k]->M); */
  /*   for(i = 0; i< U[k]->M->n; i++){ */
  /*     printf("%d ", U[k]->p[i]); */
  /*     } */
  /*   printf("\n-----------------\n"); */
  /* } */

  /* printf("-----------------\n"); */
  /* printf("L :\n"); */
  /* for(k = 0; k < diag; k++){ */
  /*   spasm_save_csr(stdout, L[k]->M); */
  /*   for(i = 0; i < L[k]->M->n; i++){ */
  /*     printf("%d ", L[k]->p[i]); */
  /*   } */
  /*   printf("\n----------------\n"); */
  /* } */

  printf("diagonales parcourues : %d sur %d\n", diag, n_blocks);
  /* calculer le rang = somme des pivots trouvés et afficher le résultat */
  rank = 0;
  for(k = 0; k < n_blocks; k++){
    rank += blocks[k].r;
  }

  printf("rang : %d\n", rank);


  /* --------- LIBÉRATION DE LA MÉMOIRE, FIN DU PROGRAMME ---------- */

  for(i = 0; i < n_blocks; i++){
    super_spasm_free(CS[i]);
    super_spasm_free(U[i]);
    free(Qinv[i]);
  }

  super_list_clear(&L_list);
  free(L);
  uptri_free(FS);
  free(CS);
  free(J);
  free(unz);
  free(un);
  free(rows);
  free(blocks);
  free(Qinv);
  free(P);
  free(U);

  return 0;
}
