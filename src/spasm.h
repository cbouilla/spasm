#ifndef _SPASM_H
#define _SPASM_H

#include <stddef.h>           // size_t
#include <inttypes.h>         // int64_t
#include <stdio.h>            // FILE
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* --- primary SpaSM routines and data structures --- */

// unfortunately we use "n" for #rows and "m" for #columns whereas the rest of the world (BLAS...)
// does the opposite... 

typedef int spasm_GFp;

typedef struct {                /* matrix in compressed-sparse row format */
	int nzmax;                    /* maximum number of entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int *p;                       /* row pointers (size n+1) */
	int *j;                       /* column indices, size nzmax */
	spasm_GFp *x;                 /* numerical values, size nzmax (optional) */
	int prime;
}      spasm;

typedef struct {                /* matrix in triplet form */
	int nzmax;                    /* maximum number of entries */
	int nz;                       /* # entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int *i;                       /* row indices, size nzmax */
	int *j;                       /* column indices (size nzmax) */
	spasm_GFp *x;                 /* numerical values, size nzmax (optional) */
	int prime;
}      spasm_triplet;


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

Compressed Row form :

p = {   0,             3,             6,        8,      10 }
i = {   0,   1,   3,   1,   2,   3,   0,   2,   1,   3 }
x = { 4.5, 3.1, 3.5, 2.9, 1.7, 0.4, 3.2, 3.0, 0.9, 1.0 }

In particular, the actual number of nnz is p[n]. Coefficients of a row need not be sorted by column index.

The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix). */

typedef struct {                /* a PLUQ factorisation */
	spasm *L;
	spasm *U;
	int *qinv;                    /* the inverse Q is stored */
	int *p;
}      spasm_lu;

typedef struct {                /* a dense LU factorization */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	int prime; 
	int *p;                       /* positions of pivots in allocated rows */
	spasm_GFp **x;                /* pointers to the rows */
}      spasm_dense_lu;


typedef struct {      /**** a Dulmage-Mendelson decomposition */
				int *p;       /* size n, row permutation */
				int *q;       /* size m, column permutation */
				int *r;       /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
				int *c;       /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
				int nb;       /* # of blocks in fine decomposition */
				int rr[5];    /* coarse row decomposition */
				int cc[5];    /* coarse column decomposition */
}      spasm_dm;


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
int spasm_nnz(const spasm * A);
void *spasm_malloc(size_t size);
void *spasm_calloc(size_t count, size_t size);
void *spasm_realloc(void *ptr, size_t size);

spasm *spasm_csr_alloc(int n, int m, int nzmax, int prime, int with_values);
void spasm_csr_realloc(spasm * A, int nzmax);
void spasm_csr_resize(spasm * A, int n, int m);
void spasm_csr_free(spasm * A);

spasm_triplet *spasm_triplet_alloc(int m, int n, int nzmax, int prime, int with_values);
void spasm_triplet_realloc(spasm_triplet * A, int nzmax);
void spasm_triplet_free(spasm_triplet * A);

spasm_dm *spasm_dm_alloc(int n, int m);
void spasm_dm_free(spasm_dm * P);

void spasm_vector_zero(spasm_GFp * x, int n);
void spasm_vector_set(spasm_GFp * x, int a, int b, spasm_GFp alpha);

spasm *spasm_identity(int n, int prime);
void spasm_human_format(int64_t n, char *target);
int spasm_get_num_threads();
int spasm_get_thread_num();

/* spasm_triplet.c */
void spasm_add_entry(spasm_triplet * T, int i, int j, spasm_GFp x);
void spasm_triplet_transpose(spasm_triplet * T);
spasm *spasm_compress(const spasm_triplet * T);

/* spasm_io.c */
spasm_triplet *spasm_load_sms(FILE * f, int prime);
spasm_triplet *spasm_load_mm(FILE * f, int prime);
spasm *spasm_load_bin(FILE * f, int prime);
void spasm_save_triplet(FILE * f, const spasm_triplet * A);
void spasm_save_csr(FILE * f, const spasm * A);
void spasm_save_pnm(const spasm * A, FILE * f, int x, int y, int mode, spasm_dm *DM);
spasm *spasm_load_gbla_old(FILE * f, int with_values);
spasm *spasm_load_gbla_new(FILE * f);

/* spasm_transpose.c */
spasm *spasm_transpose(const spasm * C, int keep_values);

/* spasm_submatrix.c */
spasm *spasm_submatrix(const spasm * A, int r_0, int r_1, int c_0, int c_1, int with_values);
spasm *sorted_spasm_submatrix(const spasm * A, int r0, int r1, int c0, int c1, int *py, int with_values);
spasm *spasm_rows_submatrix(const spasm * A, int i0, int i1, int with_values);

/* spasm_permutation.c */
void spasm_pvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
void spasm_ipvec(const int *p, const spasm_GFp * b, spasm_GFp * x, int n);
int *spasm_pinv(int const *p, int n);
spasm *spasm_permute(const spasm * A, const int *p, const int *qinv, int with_values);
int *spasm_random_permutation(int n);
void spasm_range_pvec(int *x, int a, int b, int *p);

/* spasm_GFp.c */
spasm_GFp spasm_GFp_inverse(spasm_GFp a, int prime);

/* spasm_scatter.c */
void spasm_scatter(const int *Aj, const spasm_GFp * Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime);

/* spasm_reach.c */
int spasm_dfs(int i, const spasm * G, int top, int *xi, int *pstack, int *marks, const int *pinv);
int spasm_reach(const spasm * G, const spasm * B, int k, int l, int *xi, const int *pinv);

/* spasm_gaxpy.c */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp * y);
int spasm_sparse_vector_matrix_prod(const spasm * M, const spasm_GFp * x, const int *xi, int xnz, spasm_GFp * y, int *yi);

/* spasm_triangular.c */
int spasm_is_upper_triangular(const spasm * A);
int spasm_is_lower_triangular(const spasm * A);
void spasm_dense_back_solve(const spasm * L, spasm_GFp * b, spasm_GFp * x, const int *p);
int spasm_dense_forward_solve(const spasm * U, spasm_GFp * b, spasm_GFp * x, const int *q);
int spasm_sparse_backward_solve(const spasm * L, const spasm * B, int k, int *xi, spasm_GFp * x, const int *pinv, int r_bound);
int spasm_sparse_forward_solve(const spasm * U, const spasm * B, int k, int *xi, spasm_GFp * x, const int *pinv);

/* spasm_lu.c */
spasm_lu *spasm_PLUQ(const spasm * A, const int *row_permutation, int keep_L);
spasm_lu *spasm_LU(const spasm * A, const int *row_permutation, int keep_L);
void spasm_free_LU(spasm_lu * X);
void spasm_eliminate_sparse_pivots(const spasm * A, const int npiv, const int *p, spasm_GFp *x);
int spasm_early_abort(const spasm * A, const int *p, int k, const spasm * U, int nu);

/* spasm_schur.c */
void spasm_make_pivots_unitary(spasm *A, const int *p, const int npiv);
spasm *spasm_schur(const spasm * A, const int *p, int npiv, const spasm *U, const int *qinv, double est_density, int keep_L, int *p_out);
double spasm_schur_estimate_density(const spasm * A, const int *p, int n, const spasm *U, const int *qinv, int R);
int spasm_schur_dense(const spasm *A, const int *p, int k, const spasm *U, const int *qinv, double *S, int *q);
void spasm_schur_dense_randomized(const spasm *A, const int *p, int n, const spasm *U, const int *qinv, double *S, int *q, int N, int w);

/* spasm_dense_lu.c */
spasm_dense_lu *spasm_dense_LU_alloc(int m, int prime);
void spasm_dense_LU_free(spasm_dense_lu * A);
int spasm_dense_LU_process(spasm_dense_lu *A, spasm_GFp *y);

/* spasm_solutions.c */
int spasm_PLUQ_solve(spasm * A, const spasm_GFp * b, spasm_GFp * x);
int spasm_LU_solve(spasm * A, const spasm_GFp * b, spasm_GFp * x);

/* spasm_pivots.c */
int spasm_find_pivots(spasm * A, int *p, int *qinv);
spasm * spasm_permute_pivots(const spasm *A, const int *p, int *qinv, int npiv);

/* spasm_matching.c */
int spasm_maximum_matching(const spasm * A, int *jmatch, int *imatch);
int *spasm_permute_row_matching(int n, const int *jmatch, const int *p, const int *qinv);
int *spasm_permute_column_matching(int m, const int *imatch, const int *pinv, const int *q);
int *spasm_submatching(const int *match, int a, int b, int c, int d);
int spasm_structural_rank(const spasm * A);

/* spasm_dm.c */
spasm_dm *spasm_dulmage_mendelsohn(const spasm * A);

/* spasm_scc.c */
spasm_dm *spasm_strongly_connected_components(const spasm * A);

/* spasm_cc.c */
spasm_dm *spasm_connected_components(const spasm * A, spasm * given_At);

/* spasm_kernel.c */
spasm *spasm_kernel(const spasm * A, const int *column_permutation);

/* spasm_ffpack.cpp */
void spasm_ffpack_setzero(int prime, int n, int m, double *A, int ldA);
int spasm_ffpack_echelonize(int prime, int n, int m, double *A, int ldA, size_t *Q);

/* spasm_echelonize */
struct echelonize_opts {
	/* TODO: also put here parameters of the pivot search algorithms */
	int enable_tall_and_skinny;
	int enable_dense;
	int enable_GPLU;

	/* Parameters of the "root" echelonization procedure itself */
	double min_pivot_proportion;    /* minimum number of pivots found to keep going; < 0 = keep going */
	int max_round;                  /* maximum number of rounds; < 0 = keep going */
 
 	/* Parameters that determine the choice of a finalization strategy */
 	double sparsity_threshold;      /* denser than this --> dense method; < 0 = keep going */

	/* options of dense methods */
	int dense_block_size;           /* #rows processed in each batch; determine memory consumption */
	double low_rank_ratio;          /* if k rows have rank less than k * low_rank_ratio --> "tall-and-skinny"; <0 = don't */
	double tall_and_skinny_ratio;   /* aspect ratio (#rows / #cols) higher than this --> "tall-and-skinny"; <0 = don't */
};

spasm* spasm_echelonize(spasm *A, int *Uqinv, struct echelonize_opts *opts);

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

static inline int spasm_row_weight(const spasm * A, int i) {
	int *Ap = A->p;
	return Ap[i + 1] - Ap[i];
}
#endif
