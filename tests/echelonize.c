#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <stdlib.h>

#include "spasm.h"

/* so far, this program checks that the U matrix is in REF, and that rowspan(A) is included in rowspan(U) */

void echelon_form_check(const spasm *U, int *qinv)
{
	int m = U->m;
	int *Up = U->p;
	int *Uj = U->j;
	spasm_GFp *Ux = U->x;
	
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


void deterministic_inclusion_test(const spasm *A, const spasm *U, const int *qinv)
{
	printf("# Checking that rowspan(A) is included in rowspan(U) [deterministic]...\n");
	int done = 0;
	int n = A->n;
	int m = A->m;
	#pragma omp parallel
	{
		int tid = spasm_get_thread_num();
		int *xj = spasm_malloc(3*m * sizeof(int));
	  	spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
	  	for (int j = 0; j < 3*m; j++)
	  		xj[j] = 0;
	  
	  	#pragma omp for schedule(dynamic, 10)
	  	for (int i = 0; i < n; i++) {
	  		int top = spasm_sparse_forward_solve(U, A, i, xj, x, qinv);
	  		
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


void probabilistic_inclusion_test(spasm *A, spasm *U, int n_iterations)
{
	fprintf(stderr, "---> Checking that rowspan(A) is included in rowspan(U) [probabilistic, %d iterations]...\n", n_iterations);
	
	int prime = A->prime;
	int n = A->n;
	int m = A->m;
	int r = U->n;
	int *Ap = A->p;
	int *Aj = A->j;
	spasm_GFp *Ax = A->x;
	int *Up = U->p;
	int *Uj = U->j;
	spasm_GFp *Ux = U->x;

	int done = 0;
	#pragma omp parallel
	{
		int tid = spasm_get_thread_num();
		spasm_GFp *x = spasm_malloc(m * sizeof(spasm_GFp));
	
		#pragma omp for schedule(dynamic, 1)
		for (int k = 0; k < n_iterations; k++) {
			/* x <--- random linear combination of the rows of A */
			for (int j = 0; j < m; j++)
				x[j] = 0;
			for (int i = 0; i < n; i++)
				spasm_scatter(Aj, Ax, Ap[i], Ap[i + 1], rand() % prime, x, prime);
			/* eliminate everything in x */
			for (int i = 0; i < r; i++) {
				int j = Uj[Up[i]];
				if (x[j] != 0)
					spasm_scatter(Uj, Ux, Up[i], Up[i + 1], prime - x[j], x, prime);
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
	int prime = 42013;
	int deterministic = 0;

	spasm_triplet *T = spasm_load_sms(stdin, 42013);
	spasm *A = spasm_compress(T);
	spasm_triplet_free(T);

	int m = A->m;
	int *qinv = spasm_malloc(m * sizeof(int));
	spasm *U = spasm_echelonize(A, qinv, NULL);   /* NULL = default options */
	
	assert(A->m == U->m);
	assert(U->n <= A->n);
	assert(U->n <= U->m);

	echelon_form_check(U, qinv);
	//if (deterministic)
	deterministic_inclusion_test(A, U, qinv);
	// else
	// 	probabilistic_inclusion_test(A, U, 100);
	
	// probabilistic_rank_check(A, U, qinv, 100);

	spasm_csr_free(A);
	spasm_csr_free(U);
	free(qinv);
	exit(EXIT_SUCCESS);
}
