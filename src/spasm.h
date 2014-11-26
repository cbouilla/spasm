/* indent -nfbs -i4 -nip -npsl -di0 -nut .... */

#ifndef _SPASM_H
#define _SPASM_H

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

In partiular, the actual number of nnz is p[n].

The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix).
*/

/* spasm_util.c */
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


/* spasm_triplet.c */
void spasm_add_entry(spasm_triplet *T, int i, int j, spasm_GFp x);
spasm * spasm_compress(const spasm_triplet *T);

/* spasm_io.c */
spasm_triplet * spasm_load_triplet(FILE *f, int prime);
void spasm_save_triplet(FILE *f, const spasm_triplet *A);
void spasm_save_csr(FILE *f, const spasm *A);


/* spasm_permutation.c */
#define SPASM_IDENTITY_PERMUTATION NULL

/*
void cs_pvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
void cs_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
int *spasm_pinv(int const *p, int n);
spasm *spasm_permute(const spasm * A, const int *pinv, const int *q, int values);
*/

/* spasm_GFp.c */
spasm_GFp spasm_GFp_inverse(spasm_GFp a, int prime);

/* spasm_scatter.c */
void spasm_scatter(const int *Aj, const spasm_GFp *Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime);

/* spasm_reach.c */
int spasm_dfs(int j, spasm * G, int top, int *xi, int *pstack, const int *pinv);
int spasm_reach(spasm * G, const spasm * B, int k, int *xi, const int *pinv);

/* spasm_gaxpy.c */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp *y);

/* spasm_triangular.c */
void spasm_triangular_solve(const spasm * G, spasm_GFp * x, int lo);
int spasm_sparse_triangular_solve(spasm * G, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv, int lo);



static inline int spasm_max(int a, int b) {
  return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b) {
  return (a < b) ? a : b;
}

#endif
