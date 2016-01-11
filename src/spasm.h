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

typedef struct  {
  spasm *L;
  spasm *U;
  int   *qinv;
  int   *p;
} spasm_lu;

typedef struct {
  int *p;    /* size m, row permutation */
  int *q;    /* size n, column permutation */
  int nr;    /* # of row blocks */
  int nc;    /* # of row blocks */
  int *rr;   /* row decomposition (i-th block between rr[i] and rr[i+1]) */
  int *cc;   /* column decomposition */
} spasm_partition;

typedef struct {
  spasm_partition *CC;
  spasm_partition **SCC;
} spasm_cc;


typedef struct {
  spasm_partition *DM;
  spasm_cc *H, *S, *V;
} spasm_dm;


typedef struct system {
  spasm *M; // Matrix
  spasm *B; // Right-hand side (NULL by deflaut, must be updated at each step)
  int Bnz; // number of non-zero entries in B (initialized at 0);
  int *p; // Permution over the matrix rows
  int rect; // size of the rectangular part of M
  int left; // difference of index between the current system and the system corresponding to the left part of the solution.
  int diag; // diagonal number.
  struct system *next;
}spasm_system;


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
#define SPASM_DISCARD_L 0
#define SPASM_SUCCESS 0
#define SPASM_NO_SOLUTION 1


/* spasm_util.c */
double spasm_wtime();
int spasm_nnz(const spasm *A);
void * spasm_malloc(size_t size);
void * spasm_calloc(size_t count, size_t size);
void * spasm_realloc(void *ptr, size_t size);
spasm * spasm_csr_alloc(int m, int n, int nzmax, int prime, int with_values);
void spasm_csr_realloc(spasm *A, int nzmax);
void spasm_csr_free(spasm *A);
spasm_triplet *spasm_triplet_alloc(int m, int n, int nzmax, int prime, int with_values);
void spasm_triplet_realloc(spasm_triplet *A, int nzmax);
void spasm_triplet_free(spasm_triplet *A);
void spasm_csr_resize(spasm *A, int n, int m);
spasm_partition * spasm_partition_alloc(int n, int m, int nr, int nc);
void spasm_partition_free(spasm_partition *P);
void spasm_partition_tighten(spasm_partition *P);
void spasm_vector_zero(spasm_GFp *x, int n);
void spasm_vector_set(spasm_GFp *x, int a, int b, spasm_GFp alpha);
spasm * spasm_identity(int n, int prime);


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
void spasm_save_ppm(FILE *f, const spasm *A, const spasm_dm *X);
void spasm_save_permutation(FILE *f, const int *p, int n);
int * spasm_load_permutation(FILE *f, int n);

/* spasm_transpose.c */
spasm *spasm_transpose(const spasm *C, int keep_values);

/* spasm_submatrix.c */
spasm * spasm_submatrix(const spasm *A, int r_0, int r_1, int c_0, int c_1, int with_values);
spasm * sorted_spasm_submatrix(const spasm *A, int r0, int r1, int c0, int c1, int *py, int with_values);

/* spasm_permutation.c */
void spasm_pvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
void spasm_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
int *spasm_pinv(int const *p, int n);
spasm *spasm_permute(const spasm * A, const int *pinv, const int *q, int values);
int * spasm_random_permutation(int n);
void spasm_range_pvec(int *x, int a, int b, int *p);

/* spasm_GFp.c */
spasm_GFp spasm_GFp_inverse(spasm_GFp a, int prime);

/* spasm_scatter.c */
void spasm_scatter(const int *Aj, const spasm_GFp *Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime);

/* spasm_reach.c */
int spasm_dfs(int i, const spasm * G, int top, int *xi, int *pstack, int *marks, const int *pinv);
int spasm_reach(const spasm * G, const spasm * B, int k, int l, int *xi, const int *pinv);

/* spasm_gaxpy.c */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp *y);

/* spasm_triangular.c */
int spasm_is_upper_triangular(const spasm *A);
int spasm_is_lower_triangular(const spasm *A);

void spasm_dense_back_solve(const spasm * L, spasm_GFp *b, spasm_GFp * x, const int *p);
int spasm_dense_forward_solve(const spasm * U, spasm_GFp *b, spasm_GFp * x, const int *q);
int spasm_sparse_backward_solve(const spasm * L, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv, int r_bound);
int spasm_sparse_forward_solve(const spasm * U, const spasm *B, int k, int *xi, spasm_GFp *x, const int *pinv);

/* spasm_lu.c */
spasm_lu *spasm_PLUQ(const spasm * A, const int *row_permutation, int keep_L);
spasm_lu *spasm_LU(const spasm * A, const int *row_permutation, int keep_L);
void spasm_free_LU(spasm_lu *X);
int spasm_find_pivot(int *xi, spasm_GFp *x, int top, spasm *U, spasm *L, int *unz_ptr, int *lnz_ptr, int i, int *deff_ptr, int *qinv, int *p, int n);

/* spasm_solutions.c */
int spasm_PLUQ_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x);
int spasm_LU_solve(const spasm *A, const spasm_GFp *b, spasm_GFp *x);

/* spasm_sort.c */
int * spasm_row_sort (const spasm *A);
int * spasm_cheap_pivots(const spasm *A);

/* spasm_csr.c */
void spasm_row_entries_sort(spasm *M, int with_value);

/* spasm_concatenate.c */
spasm * spasm_row_concatenation(spasm *A, spasm *B, int with_values);

/* spasm_matching.c */
int spasm_maximum_matching(const spasm *A, int *jmatch, int *imatch);
int * spasm_permute_row_matching(int n, const int *jmatch, const int *p, const int *qinv);
int * spasm_permute_column_matching(int m, const int *imatch, const int *pinv, const int *q);
int * spasm_submatching(const int *match, int a, int b, int c, int d);
int spasm_structural_rank(const spasm *A);

/* spasm_dm.c */
spasm_dm * spasm_dulmage_mendelsohn(const spasm *A);

/* spasm_cc.c */
spasm_partition * spasm_connected_components(const spasm *A, const spasm *A_t, const int *jmatch);

/* spasm_scc.c */
spasm_partition * spasm_strongly_connected_components(const spasm *A);

/* spasm_kernel.c */
spasm * spasm_kernel(const spasm *A, const int * column_permutation);

/* spasm_inverse.c */
int spasm_sparse_vector_matrix_prod(const spasm *M, const spasm_GFp *x, const int *xi, int xnz, spasm_GFp *y, int *yi);
int spasm_inverse_and_product(const spasm *L, const spasm *M, int k, spasm_GFp *y, int *yi, const int *pinv);
int spasm_solve_and_product(const spasm *L, const spasm *M, const spasm *A, int k, spasm_GFp *y, int *yi, const int *pinv);

/* spasm_lazy.c */
int spasm_add_vectors(spasm_GFp *u, int *ui, int unz, spasm_GFp *v, int *vi, int vnz, int size);
void spasm_system_right_hand_init(spasm_system *L, int nnz);
void spasm_lazy_right_hand_update(int d, int k, spasm_system *L0, spasm_system *L1, int bound, int *xi, int *x, int top, int stop, const int **p);
void spasm_new_lazy_permutation(int bound, const int *Lperm, int *p_new, int vec_size);
int spasm_next_left_system(spasm_system **L, int k, int d);
void spasm_lazy_system(spasm_system **L, int k, int l, int d, const int **p);
spasm_system * spasm_system_update(spasm_system *L, spasm *M, int *p, int rect, int left, int diag);
spasm_system * spasm_system_clear(spasm_system *L);
int spasm_lazy_computation(int d, int k, int i, spasm_system **S, spasm_GFp *u, int *ui, int usize, const spasm **A, const int **p);

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

