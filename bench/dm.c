#include <assert.h>
#include <stdio.h>
#include "spasm.h"

/* bad models :
332
175
205
457
*/



void show(spasm_cc *Y) {
  int i, j;

  for(i = 0; i < Y->CC->nr; i++) {
    int C_n = Y->CC->rr[i + 1] - Y->CC->rr[i];
    int C_m = Y->CC->cc[i + 1] - Y->CC->cc[i];
    printf("   *) Connected component (%d x %d) --- (%d, %d) to (%d, %d)\n", C_n, C_m,
	   Y->CC->rr[i], Y->CC->cc[i], Y->CC->rr[i + 1], Y->CC->cc[i + 1]);
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	int SCC_n = Y->SCC[i]->rr[j + 1] - Y->SCC[i]->rr[j];
	int SCC_m = Y->SCC[i]->cc[j + 1] - Y->SCC[i]->cc[j];
	printf("       *) SCC (%d x %d) --- (%d, %d) to (%d, %d)\n", SCC_n, SCC_m,
	       Y->SCC[i]->rr[j], Y->SCC[i]->cc[j], Y->SCC[i]->rr[j + 1], Y->SCC[i]->cc[j + 1]);
      }
    }
  }
}

int largest_diagonal_block(spasm_cc *Y) {
  int i, j, r = -1;

  for(i = 0; i < Y->CC->nr; i++) {
    if (Y->SCC[i] != NULL) {
      for(j = 0; j < Y->SCC[i]->nr; j++) {
	int SCC_n = Y->SCC[i]->rr[j + 1] - Y->SCC[i]->rr[j];
	r = spasm_max(r, SCC_n);
      }
    }
  }
  return r;
}


int main() {
    spasm_triplet *T;
    spasm *A, *B;
    spasm_dm *x;
    int n, m, i, j, *qinv;

    T = spasm_load_sms(stdin, -1);
    A = spasm_compress(T);
    spasm_triplet_free(T);

    n = A->n;
    m = A->m;

    x = spasm_dulmage_mendelsohn(A);
    i = -1;
    if (x->H != NULL) {
      i = spasm_max(i, largest_diagonal_block(x->H));
    }
    if (x->S != NULL) {
      i = spasm_max(i, largest_diagonal_block(x->S));
    }
    if (x->V != NULL) {
      i = spasm_max(i, largest_diagonal_block(x->V));
    }

    //    printf("A : %d x %d with %d nnz ----> largest block of size %.1f %%\n", n, m, spasm_nnz(A), 100.0 * i / spasm_min(n, m));
    printf("%5d \t %5d \t %6d \t %6d \t %.1f\n", n, m, spasm_nnz(A), i, 100.0 * i / spasm_min(n, m));

    /*    if (x->H != NULL) {
      int h_n = x->DM->rr[1] - x->DM->rr[0];
      int h_m = x->DM->cc[2] - x->DM->cc[0];
      printf("*) H (%d x %d) : \n", h_n, h_m);
      show(x->H);
    }
    if (x->S != NULL) {
      int s_n = x->DM->rr[2] - x->DM->rr[1];
      int s_m = x->DM->cc[3] - x->DM->cc[2];
      printf("*) S (%d x %d) : \n", s_n, s_m);
      show(x->S);
    }

    if (x->V != NULL) {
      int v_n = x->DM->rr[4] - x->DM->rr[2];
      int v_m = x->DM->cc[4] - x->DM->cc[3];
      printf("*) V (%d x %d) : \n", v_n, v_m);
      show(x->V);
      }*/


    /*
qinv = spasm_pinv(x->DM->q, m);
    B = spasm_permute(A, x->DM->p, qinv, SPASM_IGNORE_VALUES);
    free(qinv);

    
      FILE *f = fopen("permuted.sms", "w");
    spasm_save_csr(f, B);
    fclose(f);
    
    
        FILE *f = fopen("plop.ppm", "w");
    spasm_save_ppm(f, B, x);
    fclose(f);

    spasm_csr_free(B);
    */    
    spasm_csr_free(A);
    return 0;
}
