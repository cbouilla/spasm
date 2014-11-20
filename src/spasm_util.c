#include <assert.h>
#include "spasm.h"

int spasm_is_csc(const spasm *A) {
  assert (A != NULL);
  return (A->nz < 0);
}

int spasm_is_triplet(const spasm *A) {
  assert (A != NULL);
  return (A->nz >= 0);
}

int spasm_nnz(const spasm *A) {
  assert (A != NULL);
  if (spasm_is_csc(A)) {
    return A->p[A->n];
  } else {
    return (A->nz);
  }
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
  if (x == NULL) {
    perror("realloc failed");
    exit(1);
  }
  return x;
}


/* allocate a sparse matrix (triplet form or compressed-column form) */
spasm *spasm_spalloc(int m, int n, int nzmax, int prime, int values, int triplet) {
  spasm *A;

  A = spasm_malloc(sizeof(spasm)); /* allocate the cs struct */
    A->m = m;                        /* define dimensions and nzmax */
    A->n = n;
    A->nzmax = nzmax;
    A->prime = prime;
    A->nz = triplet ? 0 : -1;   /* allocate triplet or comp.col */
    A->p = spasm_malloc((triplet ? nzmax : n + 1) * sizeof(int));
    A->i = spasm_malloc(nzmax * sizeof(int));
    A->x = values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL;
    return A;
}

/* change the max # of entries in a sparse matrix.
 * If nzmax < 0, then the matrix is trimmed to its current nnz.
 */
void spasm_sprealloc(spasm *A, int nzmax) {
    assert (A != NULL);
    if (nzmax < 0) {
      nzmax = spasm_nnz(A);
    }
    A->i = spasm_realloc(A->i, nzmax * sizeof(int));
    if (spasm_is_triplet(A)) {
      A->p = realloc(A->p, nzmax * sizeof(int));
    }
    if (A->x != NULL) {
      A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
    }
    A->nzmax = nzmax;
}

/* free a sparse matrix */
void spasm_spfree(spasm *A) {
  assert(A != NULL);
  free(A->p);
  free(A->i);
  free(A->x); // trick : free does nothing on NULL pointer
  free(A);
}
