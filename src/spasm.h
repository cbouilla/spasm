/* indent -nfbs -i4 -nip -npsl -di0 -nut .... */
#ifndef _SPASM_H
#define _SPASM_H

#define SPASM_TIMING
//#define SPASM_COL_WEIGHT_PIVOT_SELECTION

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

/* --- primary SpaSM routines and data structures --- */

typedef int spasm_GFp;

typedef struct {   /* matrix in compressed-column or triplet form */
    int nzmax;     /* maximum number of entries */
    int n;         /* number of rows */
    int m;         /* number of columns */
    int *p;        /* row pointers (size n+1) */
    int *j;        /* column indices, size nzmax */
    spasm_GFp *x;  /* numerical values, size nzmax (optional) */
    int prime;
} spasm;

typedef struct  {   /* matrix in compressed-column or triplet form */
    int nzmax;      /* maximum number of entries */
    int nz;         /* # entries */
    int n;          /* number of rows */
    int m;          /* number of columns */
    int *i;         /* row indices, size nzmax */
    int *j;         /* column indices (size nzmax) */
    spasm_GFp *x;   /* numerical values, size nzmax (optional) */
    int prime;
} spasm_triplet;


/* example (this is Matrix/t1)

    [ 4.5  0.0  3.2  0.0 ]
    [ 3.1  2.9  0.0  0.9 ]
A = [ 0.0  1.7  3.0  0.0 ]
    [ 3.5  0.4  0.0  1.0 ]

Triplet form (nz != -1) :

i = {   2,   1,   3,   0,   1,   3,   3,   1,   0,   2 }
j = {   2,   0,   3,   2,   1,   0,   1,   3,   0,   1 }
x = { 3.0, 3.1, 1.0, 3.2, 2.9, 3.5, 0.4, 0.9, 4.5, 1.7 }

the coefficients may appear in any order.

Compressed Column form :

p = {   0,             3,             6,        8,      10 }
i = {   0,   1,   3,   1,   2,   3,   0,   2,   1,   3 }
x = { 4.5, 3.1, 3.5, 2.9, 1.7, 0.4, 3.2, 3.0, 0.9, 1.0 }

In particular, the actual number of nnz is p[n].

The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix).
*/

#define SPASM_IDENTITY_PERMUTATION NULL
#define SPASM_IGNORE NULL
#define SPASM_IGNORE_VALUES 0
#define SPASM_WITH_NUMERICAL_VALUES 1
#define SPASM_KEEP_L 1
#define SPASM_SUCCESS 0
#define SPASM_NO_SOLUTION 1


/* spasm_util.c */
double spasm_wtime();
int spasm_nnz(const spasm *A);

void * spasm_malloc(size_t size);
void * spasm_calloc(size_t count, size_t size);
void * spasm_realloc(void *ptr, size_t size);

spasm *spasm_csr_alloc(int m, int n, int nzmax, int prime, int with_values);
void spasm_csr_realloc(spasm *A, int nzmax);
void spasm_csr_free(spasm *A);

spasm_triplet *spasm_triplet_alloc(int m, int n, int nzmax, int prime, int with_values);
void spasm_triplet_realloc(spasm_triplet *A, int nzmax);
void spasm_triplet_free(spasm_triplet *A);

void spasm_csr_resize(spasm *A, int n, int m);

/* spasm_triplet.c */
void spasm_add_entry(spasm_triplet *T, int i, int j, spasm_GFp x);
void spasm_triplet_transpose(spasm_triplet *T);
spasm * spasm_compress(const spasm_triplet *T);

/* spasm_io.c */
spasm_triplet * spasm_load_sms(FILE *f, int prime);
void spasm_save_triplet(FILE *f, const spasm_triplet *A);
void spasm_save_csr(FILE *f, const spasm *A);
void spasm_save_pbm(FILE *f, const spasm *A);
void spasm_save_pgm(FILE *f, int x, int y, const spasm *A);
void spasm_save_ppm(FILE *f, int x, int y, const spasm *A, const spasm_partition *P);

/* spasm_transpose.c */
spasm *spasm_transpose(const spasm *C, int keep_values);

/* spasm_submatrix.c */
spasm * spasm_submatrix(const spasm *A, int r_0, int r_1, int c_0, int c_1, int with_values);


/* spasm_permutation.c */
void spasm_pvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
void spasm_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
int *spasm_pinv(int const *p, int n);
spasm *spasm_permute(const spasm * A, const int *pinv, const int *q, int values);
int * spasm_random_permutation(int n);
 
/* spasm_GFp.c */
spasm_GFp spasm_GFp_inverse(spasm_GFp a, int prime);

/* spasm_scatter.c */
void spasm_scatter(const int *Aj, const spasm_GFp *Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime);

/* spasm_reach.c */

int spasm_dfs(int i, const spasm * G, int top, int *xi, int *pstack, int *marks, const int *pinv);
int spasm_reach(const spasm * G, const spasm * B, int k, int *xi, const int *pinv);

/* spasm_gaxpy.c */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp *y);

/* spasm_triangular.c */
int spasm_is_upper_triangular(const spasm *A);
int spasm_is_lower_triangular(const spasm *A);

void spasm_dense_back_solve(const spasm * L, spasm_GFp *b, spasm_GFp * x, const int *p);
int spasm_dense_forward_solve(const spasm * U, spasm_GFp *b, spasm_GFp * x, const int *q);
int spasm_sparse_forward_solve(spasm * U, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv);

/* spasm_lu.c */
typedef struct  {
  spasm *L;
  spasm *U;
  int   *qinv;
  int   *p;
} spasm_lu;

spasm_lu *spasm_PLUQ(const spasm * A, const int *row_permutation);
spasm_lu *spasm_LU(const spasm * A, const int *row_permutation, int keep_L);
void spasm_free_LU(spasm_lu *X);

/* spasm_solutions.c */
int spasm_PLUQ_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x);
int spasm_LU_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x);

/* spasm_sort.c */
int * spasm_row_sort (const spasm *A);
int * spasm_cheap_pivots(const spasm *A);

/* spasm_matching.c */
int spasm_maximum_matching(const spasm *A, int *imatch, int *jmatch);
int spasm_structural_rank(const spasm *A);

/* spasm_dm.c */
typedef struct {
    int *p ;        /* size m, row permutation */
    int *q ;        /* size n, column permutation */
    int *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    int *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    int nb ;        /* # of blocks in fine dmperm decomposition */
    int rr [5] ;    /* coarse row decomposition */
    int cc [5] ;    /* coarse column decomposition */
} spasm_dm;

spasm_dm *spasm_dm_alloc(int n, int m);
spasm_dm * spasm_dulmage_mendelson(const spasm *A);
void spasm_dm_free(spasm_dm *D);

/* spasm_cc.c */
int spasm_connected_components(const spasm *A, int *rmark, int *cmark);

/* spasm_scc.c */
int spasm_strongly_connected_components(const spasm *A, int *p, int *r);


/* utilities */
static inline int spasm_max(int a, int b) {
  return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b) {
  return (a < b) ? a : b;
}

static inline void spasm_swap(int *a, int i, int j) {
  int x = a[i];
    a[i] = a[j];
    a[j] = x;
}

static inline int spasm_row_weight(const spasm *A, int i) {
  int *Ap;
  Ap = A->p;
  return Ap[i + 1] - Ap[i];
}


#ifdef SPASM_TIMING
#include "cycleclock.h"
#endif


#endif
