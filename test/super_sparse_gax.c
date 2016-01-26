#include <stdio.h>
#include <assert.h>
#include "spasm.h"

int main(int argc, char **argv) {
 spasm_triplet *T;
 spasm *C;
  super_spasm *M;
  spasm_GFp *x, *y, *Cx;
  int i, k, l, j, n, Cn, m, nz, xnz, *xi, *yi, test, *Mp, *Cj, *Cp;

  assert(argc > 1);
  test = atoi(argv[1]);

  T = spasm_load_sms(stdin, 277);
  C = spasm_compress(T);
  spasm_triplet_free(T);

  // créer la matrice super_spasm M.
  Cn = C->n;
  Cj = C->j;
  Cx = C->x;
  Cp = C->p;

  Mp = spasm_malloc(Cn * sizeof(int));

  M = spasm_malloc(sizeof(super_spasm));
  M->M = C;
  M->n = Cn + 3;
  M->p = Mp;
 
  // lignes non vides :
  for(i = 0; i < Cn/2; i++){
    Mp[i] = i+1;
  }
  for(i = Cn/2; i < Cn; i++){
    Mp[i] = i+2;
  }

  m = C->m; 
  n = M->n;
  xnz = 1;
  
  xi = malloc(xnz * sizeof(int));
  x = malloc(n * sizeof(spasm_GFp));
  y = malloc(C->m * sizeof(spasm_GFp));
  // z = malloc(m * sizeof(spasm_GFp));
  yi = malloc(C->m * sizeof(int));
 
  for(i = 0; i < m; i++) {
    y[i] = 0;
    // z[i] = 0;
    yi[i] = 0;
  }

  // Vecteur x:
 
  for(i = 0; i < n; i++) {
    x[i] = 0;
  }
 
  l = 0;
  for(i = 0; i < n; i++){
    x[i] = 1;
    xi[0] = i; // vecteur canonique.

    nz = super_sparse_gax(M, x, xi, xnz, y, yi);

    if(l == Cn || Mp[l] > i){ // ligne i vide.
      if(nz != 0){
	printf("not ok - empty row\n");
	exit(0);
      }
    }
    else{
      assert(Mp[l]==i); // vérifier que i existe.
      if(nz != Cp[l+1] - Cp[l]){
	printf("not ok - number of entries\n");
	exit(0);
      }

      for(k = Cp[l]; k < Cp[l+1]; k++){
	j = Cj[k];
	if(Cx[k] - y[j] != 0){
	  printf("not ok - values \n");
	  exit(0);
	}
      }
      l++;
    }
    spasm_vector_zero(x, n);
    spasm_vector_zero(xi, xnz);
    spasm_vector_zero(y, m);
    spasm_vector_zero(yi, m);
 }

 printf("ok %d - super sparse gax\n", test);
 spasm_csr_free(C);
 free(Mp);
 free(M);
 free(xi);
 free(x);
 free(yi);
 free(y);
 return 0;
}
