#ifndef _SPASM_H
#define _SPASM_H

#include <stddef.h>           // size_t
#include <inttypes.h>         // int64_t
#include <stdio.h>            // FILE
#include <stdbool.h>

typedef int64_t i64;
typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t i32;

#ifdef _OPENMP
#include <omp.h>
#endif

/* --- primary SpaSM routines and data structures --- */

// unfortunately we use "n" for #rows and "m" for #columns whereas the rest of the world (BLAS...)
// does the opposite... 

typedef i32 spasm_ZZp;

struct spasm_field_struct {
	i64 p;
    i64 halfp;
    i64 mhalfp;
    double dinvp;
};

typedef struct spasm_field_struct spasm_field[1];
// truc de la mort avec tableau de taille 1

typedef struct {                /* matrix in compressed-sparse row format */
	i64 nzmax;                    /* maximum number of entries */
	int n;                        /* number of rows */
	int m;                        /* number of columns */
	i64 *p;                       /* row pointers (size n+1) */
	int *j;                       /* column indices, size nzmax */
	spasm_ZZp *x;                 /* numerical values, size nzmax (optional) */
	spasm_field field;
	/*
	 * The actual number of entries is p[n]. 
	 * Coefficients of a row need not be sorted by column index.
	 * The numerical values are optional (useful for storing a sparse graph, or the pattern of a matrix).
	 */
} spasm;

typedef struct {                   /* matrix in triplet form */
	i64 nzmax;                     /* maximum number of entries */
	i64 nz;                        /* # entries */
	int n;                         /* number of rows */
	int m;                         /* number of columns */
	int *i;                        /* row indices, size nzmax */
	int *j;                        /* column indices (size nzmax) */
	spasm_ZZp *x;                  /* numerical values, size nzmax (optional) */
	spasm_field field;
} spasm_triplet;

typedef struct {                   /* a PLUQ factorisation */
	spasm *L;
	spasm *U;
	int *Uqinv;                    /* locate pivots in U (on column j, row Uqinv[j]) */
	int *Lqinv;                    /* locate pivots in L (on column j, row Lqinv[j]) */
	spasm_triplet *Ltmp;           /* for internal use during the factorization */
} spasm_lu;

typedef struct {      /**** a Dulmage-Mendelson decomposition */
				int *p;       /* size n, row permutation */
				int *q;       /* size m, column permutation */
				int *r;       /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
				int *c;       /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
				int nb;       /* # of blocks in fine decomposition */
				int rr[5];    /* coarse row decomposition */
				int cc[5];    /* coarse column decomposition */
} spasm_dm;

struct echelonize_opts {
	/* pivot search sub-algorithms */
	bool enable_greedy_pivot_search;

	/* echelonization sub-algorithms */
	bool enable_tall_and_skinny;
	bool enable_dense;
	bool enable_GPLU;

	/* Parameters of the "root" echelonization procedure itself */
	bool L;                         /* should we compute L / Lqinv in addition to U / Uqinv ? */
	double min_pivot_proportion;    /* minimum number of pivots found to keep going; < 0 = keep going */
	int max_round;                  /* maximum number of rounds; < 0 = keep going */
 
 	/* Parameters that determine the choice of a finalization strategy */
 	double sparsity_threshold;      /* denser than this --> dense method; < 0 = keep going */

	/* options of dense methods */
	int dense_block_size;           /* #rows processed in each batch; determine memory consumption */
	double low_rank_ratio;          /* if k rows have rank less than k * low_rank_ratio --> "tall-and-skinny"; <0 = don't */
	double tall_and_skinny_ratio;   /* aspect ratio (#rows / #cols) higher than this --> "tall-and-skinny"; <0 = don't */
	double low_rank_start_weight;   /* compute random linear combinations of this many rows; -1 = auto-select */

};

#define SPASM_IDENTITY_PERMUTATION NULL
#define SPASM_IGNORE NULL
#define SPASM_IGNORE_VALUES 0
#define SPASM_WITH_NUMERICAL_VALUES 1

/* spasm_ZZp.c */
void spasm_field_init(i64 p, spasm_field F);
spasm_ZZp spasm_ZZp_init(const spasm_field F, i64 x);
spasm_ZZp spasm_ZZp_add(const spasm_field F, spasm_ZZp a, spasm_ZZp b);
spasm_ZZp spasm_ZZp_sub(const spasm_field F, spasm_ZZp a, spasm_ZZp b);
spasm_ZZp spasm_ZZp_mul(const spasm_field F, spasm_ZZp a, spasm_ZZp b);
spasm_ZZp spasm_ZZp_inverse(const spasm_field F, spasm_ZZp a);
spasm_ZZp spasm_ZZp_axpy(const spasm_field F, spasm_ZZp a, spasm_ZZp x, spasm_ZZp y);


/* spasm_util.c */
double spasm_wtime();
i64 spasm_nnz(const spasm * A);
void *spasm_malloc(i64 size);
void *spasm_calloc(i64 count, i64 size);
void *spasm_realloc(void *ptr, i64 size);
spasm *spasm_csr_alloc(int n, int m, i64 nzmax, i64 prime, bool with_values);
void spasm_csr_realloc(spasm * A, i64 nzmax);
void spasm_csr_resize(spasm * A, int n, int m);
void spasm_csr_free(spasm * A);
spasm_triplet *spasm_triplet_alloc(int m, int n, i64 nzmax, i64 prime, bool with_values);
void spasm_triplet_realloc(spasm_triplet * A, i64 nzmax);
void spasm_triplet_free(spasm_triplet * A);
spasm_dm *spasm_dm_alloc(int n, int m);
void spasm_dm_free(spasm_dm * P);
void spasm_lu_free(spasm_lu *N);
void spasm_human_format(int64_t n, char *target);
int spasm_get_num_threads();
int spasm_get_thread_num();
static inline i64 spasm_get_prime(const spasm *A) { return A->field->p; }

/* spasm_triplet.c */
void spasm_add_entry(spasm_triplet * T, int i, int j, spasm_ZZp x);
void spasm_triplet_transpose(spasm_triplet * T);
spasm *spasm_compress(const spasm_triplet * T);

/* spasm_io.c */
spasm_triplet *spasm_load_sms(FILE * f, i64 prime);
spasm_triplet *spasm_load_mm(FILE * f, i64 prime);
void spasm_save_triplet(FILE * f, const spasm_triplet * A);
void spasm_save_csr(FILE * f, const spasm * A);
void spasm_save_pnm(const spasm * A, FILE * f, int x, int y, int mode, spasm_dm *DM);

/* spasm_transpose.c */
spasm *spasm_transpose(const spasm * C, int keep_values);

/* spasm_submatrix.c */
spasm *spasm_submatrix(const spasm * A, int r_0, int r_1, int c_0, int c_1, int with_values);

/* spasm_permutation.c */
void spasm_pvec(const int *p, const spasm_ZZp * b, spasm_ZZp * x, int n);
void spasm_ipvec(const int *p, const spasm_ZZp * b, spasm_ZZp * x, int n);
int *spasm_pinv(int const *p, int n);
spasm *spasm_permute(const spasm * A, const int *p, const int *qinv, int with_values);
int *spasm_random_permutation(int n);
void spasm_range_pvec(int *x, int a, int b, int *p);

/* spasm_scatter.c */
void spasm_scatter(const spasm *A, int i, spasm_ZZp beta, spasm_ZZp * x);

/* spasm_reach.c */
int spasm_dfs(int i, const spasm * G, int top, int *xi, int *pstack, int *marks, const int *pinv);
int spasm_reach(const spasm * A, const spasm * B, int k, int l, int *xj, const int *qinv);

/* spasm_spmv.c */
void spasm_xApy(const spasm_ZZp *x, const spasm *A, spasm_ZZp *y);
void spasm_Axpy(const spasm *A, const spasm_ZZp *x, spasm_ZZp *y);

/* spasm_triangular.c */
void spasm_dense_back_solve(const spasm *L, spasm_ZZp *b, spasm_ZZp *x, const int *p);
bool spasm_dense_forward_solve(const spasm * U, spasm_ZZp * b, spasm_ZZp * x, const int *q);
int spasm_sparse_triangular_solve(const spasm *U, const spasm *B, int k, int *xj, spasm_ZZp * x, const int *qinv);

/* spasm_schur.c */
spasm *spasm_schur(const spasm * A, const int *p, int npiv, const spasm *U, const int *qinv, double est_density, int keep_L, int *p_out);
double spasm_schur_estimate_density(const spasm * A, const int *p, int n, const spasm *U, const int *qinv, int R);
int spasm_schur_dense(const spasm *A, const int *p, int k, const spasm *U, const int *qinv, double *S, int *q);
void spasm_schur_dense_randomized(const spasm *A, const int *p, int n, const spasm *U, const int *qinv, double *S, int *q, int N, int w);

/* spasm_pivots.c */
int spasm_pivots_extract_structural(const spasm *A, spasm *U, int *Uqinv, int *p, struct echelonize_opts *opts);

/* spasm_matching.c */
int spasm_maximum_matching(const spasm *A, int *jmatch, int *imatch);
int *spasm_permute_row_matching(int n, const int *jmatch, const int *p, const int *qinv);
int *spasm_permute_column_matching(int m, const int *imatch, const int *pinv, const int *q);
int *spasm_submatching(const int *match, int a, int b, int c, int d);
int spasm_structural_rank(const spasm *A);

/* spasm_dm.c */
spasm_dm *spasm_dulmage_mendelsohn(const spasm *A);

/* spasm_scc.c */
spasm_dm *spasm_strongly_connected_components(const spasm *A);

/* spasm_ffpack.cpp */
int spasm_ffpack_echelonize(i64 prime, int n, int m, double *A, int ldA, size_t *qinv);

/* spasm_echelonize */
void spasm_echelonize_init_opts(struct echelonize_opts *opts);
spasm_lu* spasm_echelonize(const spasm *A, struct echelonize_opts *opts);

/* spasm_rref.c */
spasm * spasm_rref(const spasm_lu *fact, int *Rqinv);

/* spasm_kernel.c */
spasm * spasm_kernel(const spasm_lu *fact);
spasm * spasm_kernel_from_rref(const spasm *R, const int *qinv);

/* spasm_solve.c */
bool spasm_solve(const spasm_lu *fact, const spasm_ZZp *b, spasm_ZZp *x);

/* utilities */
static inline int spasm_max(int a, int b)
{
	return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b)
{
	return (a < b) ? a : b;
}

static inline int spasm_row_weight(const spasm * A, int i)
{
	i64 *Ap = A->p;
	return Ap[i + 1] - Ap[i];
}
#endif
