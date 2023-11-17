#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <err.h>

#include "spasm.h"

i64 prime = 42013;

void parse_command_line_options(int argc, char **argv)
{
        struct option longopts[] = {
                {"modulus", required_argument, NULL, 'p'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'p':
                        prime = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
}

/* so far, this program checks that the U matrix is in REF, and that rowspan(A) is included in rowspan(U) */
void echelon_form_check(const struct spasm_csr *U, int *qinv)
{
	int m = U->m;
	const i64 *Up = U->p;
	const int *Uj = U->j;
	const spasm_ZZp *Ux = U->x;
	
	printf("# Checking that U is really in echelon form...\n");
	for (int j = 0; j < m; j++)
		qinv[j] = -1;
	for (int i = 0; i < U->n; i++) {
		if (Up[i + 1] == Up[i])
			errx(1, "Row %d of U is empty!\n", i);
		int j = Uj[Up[i]];
		if (qinv[j] != -1)
			errx(1, "First entry of row %d is under another pivot at (%d, %d)!\n", i, qinv[j], j);
		if (Ux[Up[i]] != 1)
			errx(1, "Pivot at (%d, %d) is not unitary!\n", i, j);
		qinv[j] = i;
	}
}


void rref_check(const struct spasm_csr *U, int *qinv)
{
	const i64 *Up = U->p;
	const int *Uj = U->j;
	const spasm_ZZp *Ux = U->x;
	
	printf("# Checking that R is really in RREF...\n");
	for (int i = 0; i < U->n; i++) {
		assert(Up[i] < Up[i + 1]);           /* non-empty row */
		int j = Uj[Up[i]];
		assert(qinv[j] == i);                /* pivot is first row entry */
		assert(Ux[Up[i]] == 1);              /* pivot == 1 */
		
		/* row does not contain entries on pivotal columns */
		for (int px = Up[i] + 1; px < Up[i + 1]; px++) {
			int j = Uj[px];
			assert(qinv[j] < 0);
		}

		/* TODO : check triangular shape */
	}
}

void deterministic_inclusion_test(const struct spasm_csr *A, const struct spasm_csr *U, const int *qinv)
{
	printf("# Checking that rowspan(A) is included in rowspan(U) [deterministic]...\n");
	int done = 0;
	int n = A->n;
	int m = A->m;
	#pragma omp parallel
	{
		int tid = spasm_get_thread_num();
		int *xj = spasm_malloc(3*m * sizeof(int));
	  	spasm_ZZp *x = spasm_malloc(m * sizeof(spasm_ZZp));
	  	for (int j = 0; j < 3*m; j++)
	  		xj[j] = 0;
	  
	  	#pragma omp for schedule(dynamic, 10)
	  	for (int i = 0; i < n; i++) {
	  		int top = spasm_sparse_triangular_solve(U, A, i, xj, x, qinv);
	  		
			for (int px = top; px < m; px++) {
				int j = xj[px];
				if ((qinv[j] == -1) && (x[j] != 0))
					errx(1, "Row %d of A does not belong to rowspan(U)!\n", i);
			}

			#pragma omp atomic
			done++;

			if (tid == 0) {
				fprintf(stderr, "\r%d/%d rows checked", done, n);
				fflush(stderr);
			}

		}
		free(x);
		free(xj);
	}
	fprintf(stderr, "\n");
}


void probabilistic_inclusion_test(struct spasm_csr *A, struct spasm_csr *U, int n_iterations)
{
	fprintf(stderr, "---> Checking that rowspan(A) is included in rowspan(U) [probabilistic, %d iterations]...\n", n_iterations);
		int n = A->n;
	int m = A->m;
	int r = U->n;
	const i64 *Up = U->p;
	const int *Uj = U->j;

	int done = 0;
	#pragma omp parallel
	{
		int tid = spasm_get_thread_num();
		spasm_ZZp *x = spasm_malloc(m * sizeof(spasm_ZZp));
	
		#pragma omp for schedule(dynamic, 1)
		for (int k = 0; k < n_iterations; k++) {
			/* x <--- random linear combination of the rows of A */
			for (int j = 0; j < m; j++)
				x[j] = 0;
			for (int i = 0; i < n; i++)
				spasm_scatter(A, i, spasm_ZZp_init(A->field, rand()), x);
			/* eliminate everything in x */
			for (int i = 0; i < r; i++) {
				int j = Uj[Up[i]];
				if (x[j] != 0)
					spasm_scatter(U, i, -x[j], x);
			}
			for (int j = 0; j < m; j++)
				if ((x[j] != 0))
					errx(1, "rowspan(A) not included in rowspan(U)! (--deterministic gives a more precise diagnostic)\n");		
			done++;

			if (tid == 0) {
				fprintf(stderr, "\r%d/%d linear combinations checked", done, n_iterations);
				fflush(stderr);
			}
		}
		free(x);
	}
	fprintf(stderr, "\n");
}


/** given an arbitrary matrix A and an echelonized matrix U, check that rowspan(A) == rowspan(U). */
int main(int argc, char **argv)
{
	parse_command_line_options(argc, argv);
	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);

	int m = A->m;
	int *Rqinv = spasm_malloc(m * sizeof(int));
	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);
	opts.enable_tall_and_skinny = 1;
	struct spasm_lu *fact = spasm_echelonize(A, &opts);   /* NULL = default options */
	struct spasm_csr *U = fact->U;
	int *Uqinv = fact->qinv;

	assert(A->m == U->m);
	assert(U->n <= A->n);
	assert(U->n <= U->m);

	echelon_form_check(U, Uqinv);
	// spasm_save_csr(stdout, U);
	deterministic_inclusion_test(A, U, Uqinv);

	struct spasm_csr *R = spasm_rref(fact, Rqinv);
	rref_check(R, Rqinv);
	deterministic_inclusion_test(A, R, Rqinv);	

	// struct spasm_csr *M = spasm_trsm(U, Uqinv, A);
	// struct spasm_csr *K = spasm_kernel_from_rref(R, Rqinv);

	spasm_csr_free(A);
	spasm_csr_free(U);
	spasm_csr_free(R);
	free(Uqinv);
	free(Rqinv);
	exit(EXIT_SUCCESS);
}
