#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "spasm.h"


int main(int argc, char **argv) {
  spasm_triplet *T, *S;
  spasm  *L, *M, *I, *I1;
  spasm_GFp *y, *x, *z;
  int i, n, m, *yi, *xi, *zi, ny, nx, nz, k, test;
  char *L_filename, *M_filename;
  FILE *L_file, *M_file;

  assert(argc > 3);
  test = atoi(argv[1]);
  L_filename = argv[2];
  M_filename = argv[3];


  // load matrix
  L_file = fopen(L_filename, "r");
  if (L_file == NULL) {
    fprintf(stderr, "not ok %d - impossible to open %s\n", test, L_filename);
    exit(1);
  }
  T = spasm_load_sms(L_file, 32003);
  fclose(L_file);
  L = spasm_compress(T);
  spasm_triplet_free(T);


  //in general :
  M_file = fopen(M_filename, "r");
  if (M_file == NULL) {
    fprintf(stderr, "not ok %d - impossible to open %s\n", test, M_filename);
    exit(1);
  }
  S = spasm_load_sms(M_file, 32003);
  fclose(M_file);
  M = spasm_compress(S);
  spasm_triplet_free(S);

  assert(M->n == L->n); // check size
  assert(L->n == L->m); //lower-triangular

  m = M->m;
  n = L->n;

  // Identity :
  I = spasm_identity(n, 32003);

  I1 = spasm_identity(m, 32003);

  yi = spasm_malloc(n * sizeof(int));
  xi = spasm_malloc(n * sizeof(int));
  zi = spasm_malloc(m * sizeof(int));

  y = spasm_malloc(n * sizeof(spasm_GFp));
  x = spasm_malloc(n * sizeof(spasm_GFp));
  z = spasm_malloc(m * sizeof(spasm_GFp));



  // boucle for
  for (k = 0; k < n; k++) {
    spasm_vector_zero(x, n);
    spasm_vector_zero(xi, n);
    spasm_vector_zero(y, n);
    spasm_vector_zero(yi, n);
    spasm_vector_zero(z, n);
    spasm_vector_zero(zi, n);

    ny = spasm_inverse_and_product(L, I, I1, k, y, yi, SPASM_IDENTITY_PERMUTATION);
    //printf("non zero entries in y : %d\n", ny);
    /* for(i = 0; i < ny; i++) { */
    /*   printf("y[%d] = %d\n", yi[i], y[yi[i]]); */
    /* } */

    //check result 
    nx = spasm_sparse_vector_matrix_prod(L, y, yi, ny, x, xi);
    //printf("non zero entries in x : %d\n", n x); 
    /* for(i = 0; i < nx; i++) { */
    /*   printf("x[%d] = %d\n", xi[i], x[xi[i]]); */
    /* } */


    x[k] = x[k] - 1;

    for(i = 0; i < n; i++) {
      if(x[i] != 0) {
	printf("not ok k = %d :  x[%d] = %d\n", k, i, x[i]);
	exit(0);
      }
    }


    //printf("ok %d - inverse\n", test);
    spasm_csr_free(I);

    x = spasm_realloc(x, m * sizeof(spasm_GFp));
    xi = spasm_realloc(xi, m * sizeof(int));
  
    spasm_vector_zero(xi, m);
    spasm_vector_zero(x, m);
    spasm_vector_zero(z, m);
    spasm_vector_zero(zi, m); 

    nz = spasm_inverse_and_product(L, M, I1, k, z, zi, SPASM_IDENTITY_PERMUTATION);

    //check result

    nx = spasm_sparse_vector_matrix_prod(M, y, yi, ny, x, xi);

    for(i = 0; i < m; i++) {
      x[i] = x[i] - z[i];
      if (x[i] != 0) {
	printf("not ok k = %d : x[%d] = %d \n", k, i, x[i]);
	exit(0);
      }
    }

  }
  // fin boucle for

  printf("ok %d - inverse and product\n", test);


  spasm_csr_free(L);
  spasm_csr_free(M);
  spasm_csr_free(I);
  free(xi);
  free(x);
  free(yi);
  free(y);
  free(zi);
  free(z);

  return 0;
}
