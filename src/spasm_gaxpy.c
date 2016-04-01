#include <assert.h>
#include "spasm.h"

/*
 * (dense vector) * (sparse) Matrix
 */

/*
 * y <--- y + x*A
 */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp * y) {
    int i, n, prime;
    int *Ap, *Aj, *Ax;

    /* check inputs */
    assert(x != NULL);
    assert(y != NULL);
    assert(A != NULL);

    n = A->n;
    Ap = A->p;
    Aj = A->j;
    Ax = A->x;
    prime = A->prime;

    for (i = 0; i < n; i++) {
        spasm_scatter(Aj, Ax, Ap[i], Ap[i + 1], x[i], y, prime);
    }
}



/*
 * (sparse vector) * (sparse) Matrix Compute y = x * M, where x and M are
 * sparse.
 * 
 * The result is scattered in y, its pattern is given by yi. The return value nz
 * is the number of non-zero entries in y.
 */
int spasm_sparse_vector_matrix_prod(const spasm * M, const spasm_GFp * x, const int *xi, int xnz, spasm_GFp * y, int *yi) {
    int p, i, j, k, m, nz, Mnz, prime, *Mp, *Mj, *w;
    spasm_GFp *Mx;

    /* check inputs */
    Mnz = spasm_nnz(M);
    assert(x != NULL);
    assert(Mnz != 0);

    m = M->m;
    Mp = M->p;
    Mj = M->j;
    Mx = M->x;
    prime = M->prime;

    //get workspace, initializing w
        w = spasm_calloc(m, sizeof(int));

    /* primo, trouver support du résultat */
    nz = 0;
    for (k = 0; k < xnz; k++) {
        i = xi[k];

        for (p = Mp[i]; p < Mp[i + 1]; p++) {
            j = Mj[p];

            if (w[j] == 0) {
                w[j] = 1;
                yi[nz] = j;
                nz++;
            }
        }
    }

    /* deuxio : scatter y */
    for (k = 0; k < xnz; k++) {
        i = xi[k];
        spasm_scatter(Mj, Mx, Mp[i], Mp[i + 1], x[i], y, prime);
    }

    /* free workspace */
    free(w);
    return nz;
}

/*
 * (sparse vector) * (sparse) Matrix Compute y = x * M, where x and M are
 * sparse.
 * 
 * The result is scattered in y, its pattern is given by yi. 
 * 
 * THE xi ARRAY MUST BE SORTED IN ASCENDING ORDER.
 *
 * The return value nz is the number of non-zero entries in y.
 */
int super_sparse_gax(const super_spasm * M, const spasm_GFp * x, const int *xi, int xnz, spasm_GFp * y, int *yi) {
  int p, i, j, k, l, m, An, nz, Anz, prime, *Ap, *Aj, *w, *Mp, *wi;
  spasm_GFp *Ax;
 const spasm *A;

    /* check inputs */
    A = M->M;
    Anz = spasm_nnz(A);
    assert(x != NULL);
    assert(Anz != 0);

    m = A->m;
    An = A->n;
    Mp = M->p;
    Ap = A->p;
    Aj = A->j;
    Ax = A->x;
    prime = A->prime;

    /* get workspace, initializing w, wi */
    w = spasm_calloc(m, sizeof(int));
    wi = spasm_malloc(xnz * sizeof(int));

    for(i = 0; i < xnz; i++){
      wi[i] = -1;
    }

    /* primo, trouver support du résultat */
    nz = 0;
    l = 0;
    for (k = 0; k < xnz; k++) {
        i = xi[k];

	// Trouver ligne correspondante dans A si elle existe (non vide).
	while(l < An-1 && Mp[l] < i) l++;
	if(Mp[l] == i){ // si la ligne i est non vide dans M, correspond à la ligne l de A.
	
	  wi[k] = l; // La k-ième coordonnée de xi correspond à la ligne l de A.

	  for (p = Ap[l]; p < Ap[l + 1]; p++) {
            j = Aj[p];

            if (w[j] == 0) {
	      w[j] = 1;
	      yi[nz] = j;
	      nz++;
            }
	  }
	}
    }

    /* deuxio : scatter y */
    for (k = 0; k < xnz; k++) {
        l = wi[k];
	if(l == -1){
	  // la ligne n'existe pas dans A.
	  continue; 
	}
	i = xi[k];

        spasm_scatter(Aj, Ax, Ap[l], Ap[l + 1], x[i], y, prime);
    }

    /* free workspace */
    free(w);
    free(wi);
    return nz;
}



/*
 * (dense vector) * (sparse) Matrix Compute y = x * M, where x and y are DENSE.
 */
void super_sparse_gaxpy_dense(const super_spasm * super_M, const spasm_GFp * x, spasm_GFp * y) {
  int i, compressed_n, prime, *Ap, *Aj, *Mp;
  spasm_GFp *Ax;
 const spasm *A;

    /* check inputs */
    A = super_M->M;
    assert(x != NULL);
    
    compressed_n = A->n;
    
    Mp = super_M->p;
    Ap = A->p;
    Aj = A->j;
    Ax = A->x;
    prime = A->prime;

    for (i = 0; i < compressed_n; i++) {
      spasm_scatter(Aj, Ax, Ap[i], Ap[i + 1], x[Mp[i]], y, prime);
    }    
}


