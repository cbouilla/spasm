#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

int main(int argc, char **argv) 
{
	parse_command_line_options(argc, argv);
	struct spasm_triplet *T = spasm_triplet_load(stdin, prime, NULL);
	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);

	int n = A->n;
	int m = A->m;
 
	struct echelonize_opts opts;
	spasm_echelonize_init_opts(&opts);

	int *p = spasm_malloc(n * sizeof(*p));
	int *Uqinv = spasm_malloc(m * sizeof(*Uqinv));
	int *Lqinv = spasm_malloc(n * sizeof(*Lqinv));
	struct spasm_csr *U = spasm_csr_alloc(n, m, spasm_nnz(A), prime, true);
	U->n = 0;
	for (int j = 0; j < m; j++)
		Uqinv[j] = -1;
	for (int i = 0; i < n; i++)
		Lqinv[i] = -1;

	spasm_lu fact;
	fact.U = U;
	fact.Uqinv = Uqinv;
	fact.L = NULL;
	fact.Lqinv = Lqinv;
	fact.Ltmp = spasm_triplet_alloc(n, n, spasm_nnz(A), prime, true);

	/* find pivots, copy to U, update L */
	int npiv = spasm_pivots_extract_structural(A, NULL, &fact, p, &opts);

	/* dump pivots */
	spasm_ZZp *y = spasm_malloc(m * sizeof(*y));
	// for (int j = 0; j < m; j++) {
	// 	int i = Uqinv[j];
	// 	if (i >= 0) {
	// 		fprintf(stderr, "U[%d] eliminates column %d\n", i, j);
	// 		for (int j = 0; j < m; j++)
	// 			y[j] = 0;
	// 		spasm_scatter(U, i, 1, y);
	// 	
	// 		fprintf(stderr, "U[%d] == (", i);
	// 		for (int j = 0; j < m; j++)
	// 			fprintf(stderr, "%s%d", (j == 0) ? "" : ", ", y[j]);
	// 		fprintf(stderr, ")\n");
	// 	}
	// }
	// for (int k = 0; k < npiv; k++)
	// 	fprintf(stderr, "pivot on row A[%d]\n", p[k]);

	/* compute dense schur complement w.r.t said pivots */
	int Sm = m - npiv;
	int Sn = n - npiv;
	i64 prime = spasm_get_prime(A);
	spasm_datatype datatype = spasm_datatype_choose(prime);
	void *S = spasm_malloc(Sn * Sm * spasm_datatype_size(datatype));
	int *p_out = spasm_malloc(Sn * sizeof(*p_out));
	int *q = spasm_malloc(Sm * sizeof(*q));
	size_t *Sqinv = spasm_malloc(Sm * sizeof(*Sqinv));                   /* for FFPACK */
	size_t *Sp = spasm_malloc(Sn * sizeof(*Sp));                         /* for FFPACK */
	spasm_schur_dense(A, p + npiv, Sn, NULL, &fact, S, datatype, q, p_out);
	// for (int i = 0; i < fact.Ltmp->nz; i++)
	// 	printf("Ltmp : (%d, %d, %d)\n", fact.Ltmp->i[i], fact.Ltmp->j[i], fact.Ltmp->x[i]);
	struct spasm_csr *L = spasm_compress(fact.Ltmp);
	i64 *Lp = L->p;
	int *Lj = L->j;
	spasm_ZZp *Lx = L->x;
	assert(L->n == n);
	// printf("========================== L\n");
	// spasm_save_csr(stdout, L);
	// printf("========================== L\n");

	/* verify result */
	for (int k = 0; k < Sn; k++) {
		int i = p[npiv + k];
		assert(p_out[k] == i);
		
		/* display A[i] */
		// printf("**************** processing A[%d]\n", i);
		// for (int j = 0; j < m; j++)
		// 	y[j] = 0;
		// spasm_scatter(A, i, 1, y);
		// fprintf(stderr, "A[%d] == (", i);
		// for (int j = 0; j < m; j++)
		// 	fprintf(stderr, "%s%d", (j == 0) ? "" : ", ", y[j]);
		// fprintf(stderr, ")\n");
			
		/* start from dense row */
		for (int j = 0; j < m; j++)
			y[j] = 0;
		for (int l = 0; l < Sm; l++) {
			int j = q[l];
			y[j] = spasm_datatype_read(S, k*Sm + l, datatype);
		}
		
		// fprintf(stderr, "A[%d] == (", i);
		// for (int j = 0; j < m; j++)
		// 	fprintf(stderr, "%s%d", (j == 0) ? "" : ", ", y[j]);
		// fprintf(stderr, ")\n");

		/* add contribution from L */
		for (int px = Lp[i]; px < Lp[i + 1]; px++) {
			// fprintf(stderr, "Adding %d x U[%d]\n", Lx[px], Lj[px]);
			spasm_scatter(U, Lj[px], Lx[px], y);

			// fprintf(stderr, "A[%d] == (", i);
			// for (int j = 0; j < m; j++)
			// 	fprintf(stderr, "%s%d", (j == 0) ? "" : ", ", y[j]);
			// fprintf(stderr, ")\n");
		}
		
		/* this should be A[i] */
		// fprintf(stderr, "final check\n");
		// fprintf(stderr, "A[%d] == (", i);
		// for (int j = 0; j < m; j++)
		// 	fprintf(stderr, "%s%d", (j == 0) ? "" : ", ", y[j]);
		// fprintf(stderr, ")\n");

		spasm_scatter(A, i, -1, y);
		for (int j = 0; j < m; j++)
			assert(y[j] == 0);
	}
	exit(EXIT_SUCCESS);
}