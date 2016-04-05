#include <assert.h>
#include <stdio.h>
#include "spasm.h"

#ifdef SPASM_TIMING
extern int64 reach, scatter, data_shuffling;
#endif

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
  int *cheap_p = spasm_cheap_pivots(C); 
  LU = spasm_LU(C, cheap_p, SPASM_DISCARD_L); // on garde L

  // L*U = P*A
  // donc il faudrait que la permutation LU->P soit en réalité P-cheap[LU->P]

  int *P = LU->p;
  for(int i=0; i < C->n; i++) {
    P[i] = cheap_p[P[i]];
  }

  free(cheap_p);

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



/**************** Fonction main *********************/

int main() {
  // charge la matrice depuis l'entrée standard
  int prime = 42013;
  spasm_triplet * T = spasm_load_sms(stdin, prime);
  printf("A : %d x %d with %d nnz (density = %.3f %%) -- loaded modulo %d\n", T->n, T->m, T->nz, 100.0 * T->nz / (1.0 * T->n * T->m), prime);

  double start_time, end_time;
  start_time = spasm_wtime();

  spasm * A = spasm_compress(T);
  spasm_triplet_free(T);

  // calcule la décomposition de A
  spasm_dm * x = spasm_dulmage_mendelsohn(A);

  // B = A permutée sous forme triangulaire par blocs
  int * qinv = spasm_pinv(x->DM->q, A->m);
  spasm * B = spasm_permute(A, x->DM->p, qinv, SPASM_WITH_NUMERICAL_VALUES);
  free(qinv);

  spasm_csr_free(A);

  // calcule la liste des blocs diagonaux et leur LU (sans permutation des lignes)
  block_t *blocks;
  spasm_lu **LU;
  int n_blocks = block_list(B, x, &blocks, &LU);
  
  free(x); // le vieil objet DM : exit

  printf("blocs diagonaux : %d\n", n_blocks);

  /* Réassemble les petites matrices L, pour obtenir une grande matrice L */
  int n = B->n;
  int *P = spasm_calloc(n, sizeof(int));
  int *N = spasm_calloc(n, sizeof(int));
  int pivots = 0;
  int non_pivots = 0;

  for(int k = 0; k < n_blocks; k++) {

    int r = blocks[k].r;
    int a = blocks[k].i0;
    int block_n = blocks[k].i1 - blocks[k].i0;
    //printf("%d : (%d, %d) -- (%d, %d) {%d x %d} [rang %d]...\n", k, a, b, blocks[k].i1, blocks[k].j1, block_n, block_m, r);

    /* le seul truc sûr, c'est que ceci marque les colonnes qui contiennent des pivots. Mais sur quelles lignes seront-ils ? */
    for(int i = 0; i < r; i++) {
      P[pivots++] = a + LU[k]->p[i];
    }
    for(int i = r; i < block_n; i++) {
      N[non_pivots++] = a + LU[k]->p[i];
    }
  }
  
  for(int i=0; i<non_pivots; i++) {
    P[n - 1 - i] = N[i];
  }

  assert(pivots + non_pivots == n);
  
  printf("Go\n");

  /* reset timers */
  reach = 0;
  scatter = 0;
  data_shuffling = 0;

  spasm_lu *BIG = spasm_LU(B, P, SPASM_DISCARD_L);
  spasm *U = BIG->U;

  end_time = spasm_wtime();
  printf("\n");

  int r = U->n;
  int m = B->m;

  printf("LU factorisation (+ everything) took %.2f s\n", end_time - start_time);
  printf("U :  %d x %d with %d nnz (density = %.3f %%)\n", r, m, spasm_nnz(U), 100.0 * spasm_nnz(U) / (1.0*r*m - r*r/2.0));
  
#ifdef SPASM_TIMING
  printf("----------------------------------------\n");
  printf("reach   : %12" PRId64 "\n", reach);
  printf("scatter : %12" PRId64 "\n", scatter);
  printf("misc    : %12" PRId64 "\n", data_shuffling);
#endif
  printf("----------------------------------------\n");
  printf("rank of A = %d\n", U->n);
  
  return 0;
}
