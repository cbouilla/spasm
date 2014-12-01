#include <assert.h>
#include "spasm.h"

int spasm_nnz(const spasm *A) {
  assert (A != NULL);

  return A->p[ A->n ];
}

void * spasm_malloc(size_t size) {
  void *x = malloc(size);
  if (x == NULL) {
    perror("malloc failed");
    exit(1);
  }
  return x;
}

void * spasm_calloc(size_t count, size_t size) {
  void *x = calloc(count, size);
  if (x == NULL) {
    perror("calloc failed");
    exit(1);
  }
  return x;
}

void * spasm_realloc(void *ptr, size_t size) {
  void *x = realloc(ptr, size);
  if (ptr != NULL && x == NULL) {
    perror("realloc failed");
    exit(1);
  }
  return x;
}


/* allocate a sparse matrix (compressed-row form) */
spasm *spasm_csr_alloc(int n, int m, int nzmax, int prime, int with_values) {
  spasm *A;

    A = spasm_malloc(sizeof(spasm)); /* allocate the cs struct */
    A->m = m;                        /* define dimensions and nzmax */
    A->n = n;
    A->nzmax = nzmax;
    A->prime = prime;
    A->p = spasm_malloc((n + 1) * sizeof(int));
    A->j = spasm_malloc(nzmax * sizeof(int));
    A->x = (with_values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL);
    return A;
}


/* allocate a sparse matrix (triplet form) */
spasm_triplet *spasm_triplet_alloc(int n, int m, int nzmax, int prime, int with_values) {
  spasm_triplet *A;

    A = spasm_malloc(sizeof(spasm_triplet));
    A->m = m;
    A->n = n;
    A->nzmax = nzmax;
    A->prime = prime;
    A->nz = 0;
    A->i = spasm_malloc(nzmax * sizeof(int));
    A->j = spasm_malloc(nzmax * sizeof(int));
    A->x = (with_values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL);
    return A;
}



/* change the max # of entries in a sparse matrix.
 * If nzmax < 0, then the matrix is trimmed to its current nnz.
 */
void spasm_csr_realloc(spasm *A, int nzmax) {
    assert (A != NULL);

    if (nzmax < 0) {
      nzmax = spasm_nnz(A);
    }
    A->j = spasm_realloc(A->j, nzmax * sizeof(int));
    if (A->x != NULL) {
      A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
    }
    A->nzmax = nzmax;
}


/* change the max # of entries in a sparse matrix.
 * If nzmax < 0, then the matrix is trimmed to its current nnz.
 */
void spasm_triplet_realloc(spasm_triplet *A, int nzmax) {
    assert (A != NULL);

    if (nzmax < 0) {
      nzmax = A->nz;
    }
    A->i = spasm_realloc(A->i, nzmax * sizeof(int));
    A->j = spasm_realloc(A->j, nzmax * sizeof(int));
    if (A->x != NULL) {
      A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
    }
    A->nzmax = nzmax;
}


/* free a sparse matrix */
void spasm_csr_free(spasm *A) {
  assert(A != NULL);
  free(A->p);
  free(A->j);
  free(A->x); // trick : free does nothing on NULL pointer
  free(A);
}

void spasm_triplet_free(spasm_triplet *A) {
  assert(A != NULL);
  free(A->i);
  free(A->j);
  free(A->x); // trick : free does nothing on NULL pointer
  free(A);
}

void spasm_csr_resize(spasm *A, int n, int m) {
  int i, p, *Ap;
  assert(A != NULL);

  A->m = m;
  // in case of column shrink, check that no entries are left outside
  Ap = spasm_realloc(A->p, (n+1) * sizeof(int));

  if (A->n < n) {
    for(i = A->n; i < n + 1; i++) {
      Ap[i] = Ap[A->n];
    }
  }
  A->n = n;
}
