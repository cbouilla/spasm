#include <sys/time.h>
#include <assert.h>
#include "spasm.h"


double spasm_wtime() {
  struct timeval ts;

  gettimeofday(&ts, NULL);
  return (double) ts.tv_sec + ts.tv_usec / 1E6;
}


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
  if (ptr != NULL && x == NULL && size != 0) {
    perror("realloc failed");
    exit(1);
  }
  return x;
}


/* allocate a sparse matrix (compressed-row form) */
spasm *spasm_csr_alloc(int n, int m, int nzmax, int prime, int with_values) {
  spasm *A;

  //assert(prime <= 46337); /* to avoid arithmetic overflow */

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
  if (A == NULL) {
    return;
  }
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
  int i, *Ap;
  assert(A != NULL);

  A->m = m;
  // in case of column shrink, check that no entries are left outside
  A->p = spasm_realloc(A->p, (n+1) * sizeof(int));

  if (A->n < n) {
    Ap = A->p;
    for(i = A->n; i < n + 1; i++) {
      Ap[i] = Ap[A->n];
    }
  }
  A->n = n;
}


spasm_partition *spasm_partition_alloc(int n, int m, int nr, int nc) {
  spasm_partition *P;

    P = spasm_malloc(sizeof(spasm_partition));
    P->p = spasm_malloc(n  * sizeof(int));
    P->q = spasm_malloc(m  * sizeof(int));
    P->rr = spasm_malloc((nr+1) * sizeof(int));
    P->cc = spasm_malloc((nc+1) * sizeof(int));
    P->nr = 0;
    P->nr = 0;
    return P;
}

void spasm_partition_tighten(spasm_partition *P) {
  assert(P != NULL);
  P->rr = spasm_realloc(P->rr, (P->nr + 1) * sizeof(int));
  P->cc = spasm_realloc(P->cc, (P->nc + 1) * sizeof(int));
}


void spasm_partition_free(spasm_partition *P) {
  if (P == NULL) {
    return;
  }
  free(P->p);
  free(P->q);
  free(P->rr);
  free(P->cc);
  free(P);
}

void spasm_vector_zero(spasm_GFp *x, int n) {
  int i;

  for(i = 0; i < n; i++) {
    x[i] = 0;
  }
}

void spasm_vector_set(spasm_GFp *x, int a, int b, spasm_GFp alpha) {
  int i;

  for(i = a; i < b; i++) {
    x[i] = alpha;
  }
}

spasm * spasm_identity(int n, int prime) {
  spasm *I; 
  int i;

  I = spasm_csr_alloc(n, n, n, prime, 1);

  for(i=0; i<n; i++) {
    I->p[i] = i;
    I->j[i] = i;
    I->x[i] = 1;
  }
  I->p[n] = n;

  return I;
}

/*
 * duplicate a matrix.
 */
spasm * spasm_duplicate(const spasm * A) {
  spasm *B;
  int k, n, m, nzmax, prime, with_values;
  int *Ap, *Aj, *Bp, *Bj;
  spasm_GFp *Ax, *Bx;

  n = A->n;
  m = A->m;
  nzmax = A->nzmax;
  prime = A->prime;
  Ap = A->p;
  Aj = A->j;
  Ax = A->x;
  with_values = (Ax != NULL) ? 1 : 0;

  B = spasm_csr_alloc(n, m, nzmax, prime, with_values);
  Bp = B->p;
  Bj = B->j;
  Bx = B->x;

  for(k = 0; k <= n; k++) {
    Bp[k] = Ap[k];
  }
  for(k = 0; k < nzmax; k++) {
    Bj[k] = Aj[k];
    if(with_values) Bx[k] = Ax[k];
  }

  return B;

}
